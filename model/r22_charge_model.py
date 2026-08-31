#!/usr/bin/env python3
"""
r22_charge_model.py

House-specific R-22 charging model based on target-superheat error.

Core definitions
----------------
    V      = R-22 VSAT from LOW pressure
    SH     = SLT - V
    SHt    = 0.5 * (3*IDWB - 80 - ODDB)
    Gi     = SH - SHt

At the correct charge:
    Gi = 0

Environmental normalization
---------------------------
    Vn = V - a*(IDWB - I0) - b*(ODDB - O0)

Target model:
    Vn_target = R

so:
    VSAT_target = R + a*(IDWB - I0) + b*(ODDB - O0)

Path model:
    Gi = H(Vn)

with a monotonic, zero-preserving atan:
    H(x) = B * [atan((R-C)/W) - atan((x-C)/W)]

This guarantees:
    H(R) = 0
    dH/dx < 0 for B>0, W>0

Target outputs from today's IDWB and ODDB:
    SHt_target  = 0.5 * (3*IDWB - 80 - ODDB)
    VSAT_target = R + a*(IDWB-I0) + b*(ODDB-O0)
    LOW_target  = P_R22(VSAT_target)
    SLT_target  = VSAT_target + SHt_target

Endpoint fitting
----------------
Accepted completed recharge endpoints are the primary data for R,a,b.

1 endpoint:
    I0, O0 = that endpoint conditions
    R      = endpoint VSAT
    a=b=0
    environmental shift is not identifiable.

2 endpoints, or rank-deficient endpoint set:
    endpoint equations constrain the target surface, but cannot uniquely
    determine both a and b. Intermediate uncapped path points are used to
    choose a provisional a,b jointly with H(Vn).

>=3 independent endpoints:
    R,a,b are fitted directly by centered weighted linear least squares:
        Vj = R + a*(Ij-I0) + b*(Oj-O0)

Path points then fit only B,C,W with R,a,b fixed.

Requires:
    numpy
    scipy
    matplotlib
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import least_squares


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

VSAT_MIN = 32.0
VSAT_MAX = 40.0

# Initial, measurement-based plateau screen. Plateau points are retained but
# are not used as ordinary path-fit observations.
PLATEAU_MARGIN_F = 2.0

# Robust path outlier threshold.
OUTLIER_ABS_F = 3.0
OUTLIER_SIGMA_MULT = 3.0

# A completed recharge endpoint is accepted as a target observation when its
# measured SH error is sufficiently close to zero. Matches the deployed
# tool's own "good"/green-zone convention (index.html's ENDPOINT_G_TOL_F),
# which is looser than this script's original 2.0F port.
ENDPOINT_G_TOL_F = 3.0

# At least this many uncapped points are needed for an atan path fit.
MIN_PATH_POINTS = 3

HUBER_F_SCALE = 1.5

# Weak regularization on the path curve's B (vertical scale) and W
# (transition width), added as extra Huber-loss residuals alongside the
# real data residuals. Without this, a handful of path points clustered in
# a narrow Vn range can be fit equally well by a sane, moderate-scale curve
# OR by a degenerate one that races B to its upper bound while W collapses
# toward its floor (a tall, needle-thin "step" whose far tail still has
# just enough residual slope across that narrow window to match the data)
# -- numerically competitive, but a wildly wrong-looking curve outside it.
# Pulling B toward a modest reference value and discouraging W near zero
# costs a fit backed by many points almost nothing, while steering an
# underdetermined few-point fit away from that boundary. Ported from the
# same fix in index.html's fitPathWithFixedTarget/fitTargetPathAssisted.
PATH_B_PRIOR = 10.0
PATH_B_REG_WEIGHT = 0.15
PATH_INV_W_REG_WEIGHT = 0.5

# Approximate R-22 P/T anchors used in this analysis.
# Replace with the application's full R-22 P/T table in production.
_R22_PSIG = np.array([43.0, 57.5, 65.7, 68.6], dtype=float)
_R22_VSAT = np.array([19.9, 32.0, 38.0, 40.0], dtype=float)
_P_TO_T = PchipInterpolator(_R22_PSIG, _R22_VSAT, extrapolate=False)
_T_TO_P = PchipInterpolator(_R22_VSAT, _R22_PSIG, extrapolate=False)


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class Point:
    session_id: str
    session_type: str
    date: str
    order: int
    stage: str
    low_psig: float
    vsat_f: float
    slt_f: float
    idwb_f: float
    odbb_f: float

    sht_f: float = 0.0
    gi_f: float = 0.0
    vn_f: Optional[float] = None
    gi_pred_f: Optional[float] = None
    residual_f: Optional[float] = None
    status: str = "unknown"  # plateau, path, outlier, outside


@dataclass
class Endpoint:
    session_id: str
    date: str
    vsat_f: float
    low_psig: float
    slt_f: float
    idwb_f: float
    odbb_f: float
    gi_f: float
    accepted: bool
    reason: str


@dataclass
class TargetModel:
    I0_f: float
    O0_f: float
    R_f: float
    a_idwb: float
    b_oddb: float
    endpoint_count: int
    endpoint_rank: int
    state: str
    # state: one-endpoint, provisional-path-assisted, endpoint-regression


@dataclass
class PathModel:
    B: float
    C_f: float
    W_f: float
    rmse_f: float
    point_count: int
    available: bool


@dataclass
class Model:
    target: TargetModel
    path: Optional[PathModel]


# ---------------------------------------------------------------------------
# Refrigerant and superheat equations
# ---------------------------------------------------------------------------

def pressure_to_vsat_r22(psig: float) -> float:
    v = float(_P_TO_T(psig))
    if not math.isfinite(v):
        raise ValueError(
            f"LOW={psig:g} psig is outside the built-in R-22 P/T range."
        )
    return v


def vsat_to_pressure_r22(vsat_f: float) -> float:
    p = float(_T_TO_P(vsat_f))
    if not math.isfinite(p):
        raise ValueError(
            f"VSAT={vsat_f:g} F is outside the built-in R-22 P/T range."
        )
    return p


def target_superheat(idwb_f: float, odbb_f: float) -> float:
    return 0.5 * (3.0 * idwb_f - 80.0 - odbb_f)


def measured_gi(vsat_f: float, slt_f: float, idwb_f: float, odbb_f: float) -> float:
    sh = slt_f - vsat_f
    return sh - target_superheat(idwb_f, odbb_f)


# ---------------------------------------------------------------------------
# JSON loading
# ---------------------------------------------------------------------------

def _f(v) -> Optional[float]:
    if v in (None, ""):
        return None
    return float(v)


def load_points(path: Path) -> list[Point]:
    with path.open("r", encoding="utf-8") as f:
        root = json.load(f)

    points: list[Point] = []

    for site in root.get("sites", []):
        for rec in site.get("archive", []):
            base = rec.get("baseline", {})
            idwb = _f(base.get("indoorWB"))
            odbb = _f(base.get("outdoorDB"))
            if idwb is None or odbb is None:
                continue

            sid = str(rec.get("id", ""))
            typ = str(rec.get("type", ""))
            ts = base.get("timestampStr") or rec.get("closedAtStr") or ""
            date = ts.split(",")[0] if ts else ""

            samples = [("baseline", base)]
            samples.extend(
                (f"step {i}", s)
                for i, s in enumerate(rec.get("steps", []), start=1)
            )

            for order, (stage, sample) in enumerate(samples):
                low = _f(sample.get("lowPsig"))
                slt = _f(sample.get("suctionTemp"))
                if low is None or slt is None:
                    continue

                try:
                    vsat = pressure_to_vsat_r22(low)
                except ValueError:
                    continue

                sht = target_superheat(idwb, odbb)
                gi = measured_gi(vsat, slt, idwb, odbb)

                points.append(
                    Point(
                        session_id=sid,
                        session_type=typ,
                        date=date,
                        order=order,
                        stage=stage,
                        low_psig=low,
                        vsat_f=vsat,
                        slt_f=slt,
                        idwb_f=idwb,
                        odbb_f=odbb,
                        sht_f=sht,
                        gi_f=gi,
                    )
                )
    return points


# ---------------------------------------------------------------------------
# Endpoint extraction and target model
# ---------------------------------------------------------------------------

def extract_recharge_endpoints(points: list[Point]) -> list[Endpoint]:
    """
    Use the last recorded point of each recharge session as the candidate
    completed-recharge endpoint.

    It is accepted only when:
      - VSAT is inside [VSAT_MIN, VSAT_MAX]
      - |Gi| <= ENDPOINT_G_TOL_F
    """
    by_session: dict[str, list[Point]] = {}
    for p in points:
        if p.session_type == "recharge":
            by_session.setdefault(p.session_id, []).append(p)

    endpoints: list[Endpoint] = []

    for sid, pts in by_session.items():
        p = sorted(pts, key=lambda x: x.order)[-1]

        if not (VSAT_MIN <= p.vsat_f <= VSAT_MAX):
            accepted = False
            reason = "VSAT outside normal target range"
        elif abs(p.gi_f) > ENDPOINT_G_TOL_F:
            accepted = False
            reason = f"|Gi|>{ENDPOINT_G_TOL_F:g} F"
        else:
            accepted = True
            reason = "accepted target endpoint"

        endpoints.append(
            Endpoint(
                session_id=sid,
                date=p.date,
                vsat_f=p.vsat_f,
                low_psig=p.low_psig,
                slt_f=p.slt_f,
                idwb_f=p.idwb_f,
                odbb_f=p.odbb_f,
                gi_f=p.gi_f,
                accepted=accepted,
                reason=reason,
            )
        )

    return sorted(endpoints, key=lambda e: e.date)


def endpoint_pivots(endpoints: list[Endpoint]) -> tuple[float, float]:
    accepted = [e for e in endpoints if e.accepted]
    if not accepted:
        raise ValueError("No accepted recharge endpoints.")
    # Equal endpoint weights for now.
    I0 = float(np.mean([e.idwb_f for e in accepted]))
    O0 = float(np.mean([e.odbb_f for e in accepted]))
    return I0, O0


def endpoint_design(endpoints: list[Endpoint], I0: float, O0: float):
    accepted = [e for e in endpoints if e.accepted]
    X = np.array(
        [[1.0, e.idwb_f - I0, e.odbb_f - O0] for e in accepted],
        dtype=float,
    )
    y = np.array([e.vsat_f for e in accepted], dtype=float)
    return accepted, X, y


def path_candidate_points(points: list[Point]) -> list[Point]:
    """
    Initial path candidates:
      - inside normal VSAT range
      - not obviously IDWB-limited

    Recharge-sequence points are used as the initial seed for provisional
    path-assisted fitting. Checkup points are admitted only after they agree
    with that seed model. This prevents a single questionable checkup from
    steering an underdetermined two-endpoint fit.
    """
    result = []
    for p in points:
        if not (VSAT_MIN <= p.vsat_f <= VSAT_MAX):
            p.status = "outside"
            continue

        if p.slt_f >= p.idwb_f - PLATEAU_MARGIN_F:
            p.status = "plateau"
            continue

        p.status = "path"
        result.append(p)

    return result


def recharge_path_seed(candidates: list[Point]) -> list[Point]:
    return [p for p in candidates if p.session_type == "recharge"]


def H(vn, B, C, W, R):
    """
    Monotonic zero-preserving path curve.

    H(R)=0 exactly.
    For B>0 and W>0, dH/dVn < 0.
    """
    x = np.asarray(vn, dtype=float)
    return B * (
        np.arctan((R - C) / W)
        - np.arctan((x - C) / W)
    )


def normalized_vsat(vsat, idwb, odbb, target: TargetModel):
    return (
        np.asarray(vsat, dtype=float)
        - target.a_idwb * (np.asarray(idwb, dtype=float) - target.I0_f)
        - target.b_oddb * (np.asarray(odbb, dtype=float) - target.O0_f)
    )


def _robust_sigma(values: np.ndarray) -> float:
    if values.size == 0:
        return 0.0
    med = float(np.median(values))
    mad = float(np.median(np.abs(values - med)))
    return 1.4826 * mad


def _fit_path_with_fixed_target(
    candidates: list[Point],
    target: TargetModel,
) -> Optional[PathModel]:
    if len(candidates) < MIN_PATH_POINTS:
        return None

    V = np.array([p.vsat_f for p in candidates])
    I = np.array([p.idwb_f for p in candidates])
    O = np.array([p.odbb_f for p in candidates])
    y = np.array([p.gi_f for p in candidates])

    Vn = normalized_vsat(V, I, O, target)

    def residual(q):
        B, C, W = q
        rr = H(Vn, B, C, W, target.R_f) - y
        reg = np.array([PATH_B_REG_WEIGHT * (B - PATH_B_PRIOR), PATH_INV_W_REG_WEIGHT / W])
        return np.concatenate([rr, reg])

    res = least_squares(
        residual,
        x0=np.array([10.0, 37.0, 1.2]),
        bounds=([0.01, VSAT_MIN, 0.15], [100.0, VSAT_MAX, 10.0]),
        loss="huber",
        f_scale=HUBER_F_SCALE,
        max_nfev=100000,
    )

    B, C, W = res.x
    pred = H(Vn, B, C, W, target.R_f)
    raw = y - pred

    sigma = _robust_sigma(raw)
    threshold = max(OUTLIER_ABS_F, OUTLIER_SIGMA_MULT * sigma)

    keep = np.abs(raw) <= threshold

    # Refit after robust rejection when enough points remain.
    if np.sum(keep) >= MIN_PATH_POINTS and not np.all(keep):
        Vn2 = Vn[keep]
        y2 = y[keep]

        def residual2(q):
            B2, C2, W2 = q
            rr = H(Vn2, B2, C2, W2, target.R_f) - y2
            reg = np.array([PATH_B_REG_WEIGHT * (B2 - PATH_B_PRIOR), PATH_INV_W_REG_WEIGHT / W2])
            return np.concatenate([rr, reg])

        res = least_squares(
            residual2,
            res.x,
            bounds=([0.01, VSAT_MIN, 0.15], [100.0, VSAT_MAX, 10.0]),
            loss="huber",
            f_scale=HUBER_F_SCALE,
            max_nfev=100000,
        )
        B, C, W = res.x

    # Update statuses/predictions.
    errs = []
    for p in candidates:
        vn = float(normalized_vsat(p.vsat_f, p.idwb_f, p.odbb_f, target))
        predp = float(H(vn, B, C, W, target.R_f))
        err = p.gi_f - predp
        p.vn_f = vn
        p.gi_pred_f = predp
        p.residual_f = err

    all_err = np.array([p.residual_f for p in candidates], dtype=float)
    sigma = _robust_sigma(all_err)
    threshold = max(OUTLIER_ABS_F, OUTLIER_SIGMA_MULT * sigma)

    kept_err = []
    kept_count = 0
    for p in candidates:
        if abs(float(p.residual_f)) > threshold:
            p.status = "outlier"
        else:
            p.status = "path"
            kept_err.append(float(p.residual_f))
            kept_count += 1

    if kept_count < MIN_PATH_POINTS:
        return None

    rmse = float(np.sqrt(np.mean(np.square(kept_err))))
    return PathModel(
        B=float(B),
        C_f=float(C),
        W_f=float(W),
        rmse_f=rmse,
        point_count=kept_count,
        available=True,
    )


def _fit_target_path_assisted(
    endpoints: list[Endpoint],
    candidates: list[Point],
    I0: float,
    O0: float,
) -> tuple[TargetModel, Optional[PathModel]]:
    """
    Used when endpoint regression cannot uniquely determine R,a,b.

    Accepted endpoint equations receive strong residual weight. Path data then
    select a provisional solution among the otherwise underdetermined target
    surfaces.

    This is intentionally marked provisional.
    """
    accepted = [e for e in endpoints if e.accepted]

    if not accepted:
        raise ValueError("No accepted endpoints.")

    # In the underdetermined/provisional state, first establish H(Vn) from
    # recharge-sequence observations. Other checkups are validation candidates
    # and are admitted later only if consistent with that provisional curve.
    seed_candidates = recharge_path_seed(candidates)
    if len(seed_candidates) >= MIN_PATH_POINTS:
        fit_candidates = seed_candidates
    else:
        fit_candidates = candidates

    if len(accepted) == 1:
        e = accepted[0]
        target = TargetModel(
            I0_f=e.idwb_f,
            O0_f=e.odbb_f,
            R_f=e.vsat_f,
            a_idwb=0.0,
            b_oddb=0.0,
            endpoint_count=1,
            endpoint_rank=1,
            state="one-endpoint",
        )
        path = _fit_path_with_fixed_target(candidates, target)
        return target, path

    # Variables: R,a,b,B,C,W
    def residual(q):
        R, a, b, B, C, W = q
        rr = []

        # Endpoint constraints: strong but not mathematically exact because
        # analog gauge/SLT measurements have finite error.
        for e in accepted:
            vt = R + a*(e.idwb_f-I0) + b*(e.odbb_f-O0)
            rr.append(4.0 * (vt - e.vsat_f))

        # Path residuals from the recharge-sequence seed.
        for p in fit_candidates:
            vn = p.vsat_f - a*(p.idwb_f-I0) - b*(p.odbb_f-O0)
            rr.append(H(vn, B, C, W, R) - p.gi_f)

        # Weak regularization only to resolve severe underdetermination.
        rr.append(0.15 * a)
        rr.append(0.15 * b)
        # Same B/W regularization as _fit_path_with_fixed_target -- matters
        # even more here, since R,a,b,B,C,W are all being fit jointly: with
        # only two endpoints and a handful of path points, that's a more
        # underdetermined problem than fitting B,C,W alone against a fixed
        # target.
        rr.append(PATH_B_REG_WEIGHT * (B - PATH_B_PRIOR))
        rr.append(PATH_INV_W_REG_WEIGHT / W)
        return np.asarray(rr)

    R0 = float(np.mean([e.vsat_f for e in accepted]))
    q0 = np.array([R0, 0.0, 0.0, 10.0, 37.0, 1.2])

    res = least_squares(
        residual,
        q0,
        bounds=(
            [VSAT_MIN, -3.0, -1.0, 0.01, VSAT_MIN, 0.15],
            [VSAT_MAX,  3.0,  1.0, 100.0, VSAT_MAX, 10.0],
        ),
        loss="huber",
        f_scale=HUBER_F_SCALE,
        max_nfev=200000,
    )

    R, a, b, B, C, W = res.x

    _, X, _ = endpoint_design(endpoints, I0, O0)
    rank = int(np.linalg.matrix_rank(X))

    target = TargetModel(
        I0_f=I0,
        O0_f=O0,
        R_f=float(R),
        a_idwb=float(a),
        b_oddb=float(b),
        endpoint_count=len(accepted),
        endpoint_rank=rank,
        state="provisional-path-assisted",
    )

    # First calculate a robust residual scale from the recharge seed only.
    seed_errors = []
    for p in fit_candidates:
        vn = float(normalized_vsat(p.vsat_f, p.idwb_f, p.odbb_f, target))
        pred = float(H(vn, B, C, W, target.R_f))
        seed_errors.append(p.gi_f - pred)

    seed_errors = np.asarray(seed_errors, dtype=float)
    seed_sigma = _robust_sigma(seed_errors)
    seed_threshold = max(OUTLIER_ABS_F, OUTLIER_SIGMA_MULT * seed_sigma)

    # Admit non-seed checkup points only when they agree with the recharge-derived
    # provisional curve. This makes "every point can add information" conditional
    # on consistency rather than allowing one checkup to define the curve.
    admitted = []
    for p in candidates:
        vn = float(normalized_vsat(p.vsat_f, p.idwb_f, p.odbb_f, target))
        pred = float(H(vn, B, C, W, target.R_f))
        err = p.gi_f - pred
        p.vn_f = vn
        p.gi_pred_f = pred
        p.residual_f = err
        if p in fit_candidates or abs(err) <= seed_threshold:
            admitted.append(p)
        else:
            p.status = "outlier"

    candidates = admitted

    # Robustly classify admitted path points against provisional model.
    for p in candidates:
        p.vn_f = float(normalized_vsat(p.vsat_f, p.idwb_f, p.odbb_f, target))
        p.gi_pred_f = float(H(p.vn_f, B, C, W, R))
        p.residual_f = p.gi_f - p.gi_pred_f

    errs = np.array([p.residual_f for p in candidates], dtype=float)
    sigma = _robust_sigma(errs)
    threshold = max(OUTLIER_ABS_F, OUTLIER_SIGMA_MULT * sigma)

    keep = []
    for p in candidates:
        if abs(float(p.residual_f)) > threshold:
            p.status = "outlier"
        else:
            p.status = "path"
            keep.append(p)

    # One final joint refit with outliers removed.
    if len(keep) >= MIN_PATH_POINTS and len(keep) < len(candidates):
        candidates = keep

        def residual2(q):
            R2, a2, b2, B2, C2, W2 = q
            rr = []
            for e in accepted:
                vt = R2 + a2*(e.idwb_f-I0) + b2*(e.odbb_f-O0)
                rr.append(4.0 * (vt - e.vsat_f))
            for p in candidates:
                vn = p.vsat_f - a2*(p.idwb_f-I0) - b2*(p.odbb_f-O0)
                rr.append(H(vn, B2, C2, W2, R2) - p.gi_f)
            rr.append(0.15*a2)
            rr.append(0.15*b2)
            return np.asarray(rr)

        res = least_squares(
            residual2,
            res.x,
            bounds=(
                [VSAT_MIN, -3.0, -1.0, 0.01, VSAT_MIN, 0.15],
                [VSAT_MAX,  3.0,  1.0, 100.0, VSAT_MAX, 10.0],
            ),
            loss="huber",
            f_scale=HUBER_F_SCALE,
            max_nfev=200000,
        )
        R, a, b, B, C, W = res.x
        target.R_f = float(R)
        target.a_idwb = float(a)
        target.b_oddb = float(b)

    path_err = []
    for p in candidates:
        p.vn_f = float(normalized_vsat(p.vsat_f, p.idwb_f, p.odbb_f, target))
        p.gi_pred_f = float(H(p.vn_f, B, C, W, target.R_f))
        p.residual_f = p.gi_f - p.gi_pred_f
        path_err.append(p.residual_f)

    path = PathModel(
        B=float(B),
        C_f=float(C),
        W_f=float(W),
        rmse_f=float(np.sqrt(np.mean(np.square(path_err)))) if path_err else math.nan,
        point_count=len(candidates),
        available=len(candidates) >= MIN_PATH_POINTS,
    )
    return target, path


def fit_model(points: list[Point]) -> tuple[Model, list[Endpoint]]:
    endpoints = extract_recharge_endpoints(points)
    accepted = [e for e in endpoints if e.accepted]

    if not accepted:
        raise ValueError(
            "No accepted completed recharge endpoints. "
            "A target surface cannot yet be established."
        )

    candidates = path_candidate_points(points)
    I0, O0 = endpoint_pivots(endpoints)

    accepted2, X, y = endpoint_design(endpoints, I0, O0)
    rank = int(np.linalg.matrix_rank(X))

    # Full endpoint regression is possible only if rank=3.
    if len(accepted2) >= 3 and rank == 3:
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        R, a, b = beta

        target = TargetModel(
            I0_f=I0,
            O0_f=O0,
            R_f=float(R),
            a_idwb=float(a),
            b_oddb=float(b),
            endpoint_count=len(accepted2),
            endpoint_rank=rank,
            state="endpoint-regression",
        )
        path = _fit_path_with_fixed_target(candidates, target)
        return Model(target=target, path=path), endpoints

    # 1 endpoint or rank-deficient set (including the common 2-endpoint case).
    target, path = _fit_target_path_assisted(
        endpoints, candidates, I0, O0
    )
    return Model(target=target, path=path), endpoints


# ---------------------------------------------------------------------------
# Target calculations
# ---------------------------------------------------------------------------

def target_vsat(idwb_f: float, odbb_f: float, target: TargetModel) -> float:
    return (
        target.R_f
        + target.a_idwb*(idwb_f-target.I0_f)
        + target.b_oddb*(odbb_f-target.O0_f)
    )


def target_values(idwb_f: float, odbb_f: float, model: Model) -> dict:
    sht = target_superheat(idwb_f, odbb_f)
    vt = target_vsat(idwb_f, odbb_f, model.target)
    low = vsat_to_pressure_r22(vt)
    slt = vt + sht

    return {
        "IDWB_f": idwb_f,
        "ODDB_f": odbb_f,
        "SHt_f": sht,
        "VSAT_target_f": vt,
        "LOW_target_psig": low,
        "SLT_target_f": slt,
    }


def predicted_gi(
    vsat_f: float, idwb_f: float, odbb_f: float, model: Model
) -> Optional[float]:
    if model.path is None or not model.path.available:
        return None
    vn = float(normalized_vsat(vsat_f, idwb_f, odbb_f, model.target))
    return float(H(
        vn,
        model.path.B,
        model.path.C_f,
        model.path.W_f,
        model.target.R_f,
    ))


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def print_report(model: Model, endpoints: list[Endpoint], points: list[Point]) -> None:
    t = model.target

    print("TARGET MODEL")
    print(f"  state             = {t.state}")
    print(f"  accepted endpoints= {t.endpoint_count}")
    print(f"  endpoint rank     = {t.endpoint_rank}")
    print(f"  I0                = {t.I0_f:.4f} F")
    print(f"  O0                = {t.O0_f:.4f} F")
    print(f"  R                 = {t.R_f:.4f} F")
    print(f"  a_IDWB            = {t.a_idwb:.6f} F_VSAT/F_IDWB")
    print(f"  b_ODDB            = {t.b_oddb:.6f} F_VSAT/F_ODDB")

    print("\nTARGET EQUATIONS")
    print("  SHt = 0.5*(3*IDWB - 80 - ODDB)")
    print("  Vn  = VSAT - a*(IDWB-I0) - b*(ODDB-O0)")
    print("  target: Vn = R")
    print("  VSAT_t = R + a*(IDWB-I0) + b*(ODDB-O0)")
    print("  LOW_t  = P_R22(VSAT_t)")
    print("  SLT_t  = VSAT_t + SHt")

    print("\nRECHARGE ENDPOINTS")
    for e in endpoints:
        mark = "ACCEPT" if e.accepted else "REJECT"
        print(
            f"  {e.date:10s} VSAT={e.vsat_f:6.3f} "
            f"LOW={e.low_psig:5.1f} SLT={e.slt_f:5.1f} "
            f"WB={e.idwb_f:5.1f} OD={e.odbb_f:5.1f} "
            f"Gi={e.gi_f:+6.3f}  {mark}: {e.reason}"
        )

    if model.path is None or not model.path.available:
        print("\nPATH MODEL")
        print("  unavailable: fewer than 3 accepted non-plateau path points")
    else:
        p = model.path
        print("\nPATH MODEL")
        print(f"  Gi = H(Vn)")
        print(
            "  H(x) = B*[atan((R-C)/W) - atan((x-C)/W)]"
        )
        print(f"  B                 = {p.B:.6f}")
        print(f"  C                 = {p.C_f:.6f} F")
        print(f"  W                 = {p.W_f:.6f} F")
        print(f"  path points       = {p.point_count}")
        print(f"  path RMSE         = {p.rmse_f:.4f} F")

    print("\nPOINT CLASSIFICATION")
    for p in points:
        print(
            f"  {p.date:10s} {p.stage:8s} "
            f"V={p.vsat_f:5.2f} SLT={p.slt_f:5.1f} "
            f"WB={p.idwb_f:5.1f} OD={p.odbb_f:5.1f} "
            f"Gi={p.gi_f:+6.2f} {p.status}"
        )


def make_plot(model: Model, endpoints: list[Endpoint], points: list[Point], path: Path) -> None:
    fig, ax = plt.subplots(figsize=(11, 7))

    plotted = [p for p in points if VSAT_MIN <= p.vsat_f <= VSAT_MAX]

    for status, marker, label in [
        ("path", "o", "Path points"),
        ("plateau", "s", "IDWB plateau"),
        ("outlier", "x", "Outliers"),
    ]:
        g = [p for p in plotted if p.status == status]
        if g:
            ax.scatter(
                [p.vsat_f for p in g],
                [p.gi_f for p in g],
                marker=marker,
                s=75 if marker != "x" else 110,
                label=label,
                zorder=3,
            )

    for p in plotted:
        ax.annotate(
            f"{p.date} {p.stage} [{p.status}]",
            (p.vsat_f, p.gi_f),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=7.0,
        )

    ax.axhline(0.0, linestyle=":", linewidth=2, label="Target Gi=0")

    # Show path curves for each accepted endpoint's environmental conditions.
    if model.path is not None and model.path.available:
        grid = np.linspace(VSAT_MIN, VSAT_MAX, 500)

        for e in [x for x in endpoints if x.accepted]:
            vn = normalized_vsat(
                grid,
                np.full_like(grid, e.idwb_f),
                np.full_like(grid, e.odbb_f),
                model.target,
            )
            y = H(
                vn,
                model.path.B,
                model.path.C_f,
                model.path.W_f,
                model.target.R_f,
            )
            mask = (y <= 30.0) & (y >= -5.0)
            ax.plot(
                grid[mask],
                y[mask],
                linestyle="--",
                linewidth=1.8,
                label=f"Fit {e.date}: WB {e.idwb_f:g}, OD {e.odbb_f:g}",
            )

    ax.set_xlim(VSAT_MIN - 0.2, VSAT_MAX + 0.2)
    ax.set_xlabel("R-22 VSAT (F)")
    ax.set_ylabel("Gi = SH - SHtarget (F)")
    ax.set_title("R-22 target-superheat error model")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def make_normalized_plot(model: Model, endpoints: list[Endpoint], points: list[Point], path: Path) -> None:
    """
    Plot Gi = H(Vn) in the horizontally normalized coordinate.

    This reproduces the preferred diagnostic view:
      - all normalized points
      - recharge sequences connected
      - common monotonic H(Vn) curve
      - horizontal Gi=0 target
      - vertical Vn=R normalized target
    """
    if model.path is None or not model.path.available:
        return

    fig, ax = plt.subplots(figsize=(11, 7))

    plotted = [p for p in points if VSAT_MIN <= p.vsat_f <= VSAT_MAX]

    # Ensure Vn is available for every point in the modeled range, including
    # plateau points that were not used by the path fit.
    for p in plotted:
        if p.vn_f is None:
            p.vn_f = float(normalized_vsat(
                p.vsat_f, p.idwb_f, p.odbb_f, model.target
            ))

    ax.scatter(
        [p.vn_f for p in plotted],
        [p.gi_f for p in plotted],
        s=65,
        label="Measured points in normalized coordinate",
        zorder=3,
    )

    # Connect recharge sequences in acquisition order.
    recharge_sessions = {}
    for p in plotted:
        if p.session_type == "recharge":
            recharge_sessions.setdefault(p.session_id, []).append(p)

    for sid, seq in recharge_sessions.items():
        seq = sorted(seq, key=lambda x: x.order)
        date = seq[0].date if seq else sid
        ax.plot(
            [p.vn_f for p in seq],
            [p.gi_f for p in seq],
            marker="o",
            linewidth=2.0,
            label=f"{date} recharge",
            zorder=2,
        )

    # Common normalized path curve.
    vns = [p.vn_f for p in plotted if p.vn_f is not None]
    xmin = min(vns) - 0.5
    xmax = max(max(vns) + 0.5, model.target.R_f + 0.6)
    grid = np.linspace(xmin, xmax, 600)
    y = H(
        grid,
        model.path.B,
        model.path.C_f,
        model.path.W_f,
        model.target.R_f,
    )
    mask = (y <= 30.0) & (y >= -5.0)
    ax.plot(
        grid[mask], y[mask], linewidth=2.5, label="Common monotonic atan H(Vn)"
    )

    ax.axhline(0.0, linestyle=":", linewidth=2.0, label="Target Gi=0")
    ax.axvline(
        model.target.R_f,
        linestyle=":",
        linewidth=1.8,
        label=f"Common normalized target Vn={model.target.R_f:.2f}",
    )

    for p in plotted:
        ax.annotate(
            f"{p.date} {p.stage} [{p.status}]",
            (p.vn_f, p.gi_f),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=7.0,
        )

    ax.set_xlabel("Normalized VSAT, Vn (F)")
    ax.set_ylabel("Gi = SH - SHtarget (F)")
    ax.set_title("Horizontal normalization: Gi vs Vn")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("json_file", type=Path)
    parser.add_argument("--plot", type=Path, default=Path("r22_fit.png"))
    parser.add_argument(
        "--normalized-plot",
        type=Path,
        help="Output path for the Gi=H(Vn) horizontal-normalization plot. "
             "Defaults to <plot-stem>_Gi_vs_Vn.png.",
    )
    parser.add_argument("--idwb", type=float)
    parser.add_argument("--oddb", type=float)
    args = parser.parse_args()

    if (args.idwb is None) != (args.oddb is None):
        raise SystemExit("--idwb and --oddb must be supplied together.")

    points = load_points(args.json_file)
    if not points:
        raise SystemExit("No usable points found.")

    model, endpoints = fit_model(points)
    print_report(model, endpoints, points)

    make_plot(model, endpoints, points, args.plot)
    print(f"\nPlot written to: {args.plot}")

    normalized_plot = args.normalized_plot
    if normalized_plot is None:
        normalized_plot = args.plot.with_name(args.plot.stem + "_Gi_vs_Vn" + args.plot.suffix)
    make_normalized_plot(model, endpoints, points, normalized_plot)
    if model.path is not None and model.path.available:
        print(f"Normalized Gi=H(Vn) plot written to: {normalized_plot}")

    if args.idwb is not None:
        v = target_values(args.idwb, args.oddb, model)
        print("\nTARGET FOR CURRENT CONDITIONS")
        print(f"  IDWB        = {v['IDWB_f']:.3f} F")
        print(f"  ODDB        = {v['ODDB_f']:.3f} F")
        print(f"  SH target   = {v['SHt_f']:.3f} F")
        print(f"  VSAT target = {v['VSAT_target_f']:.3f} F")
        print(f"  LOW target  = {v['LOW_target_psig']:.3f} psig")
        print(f"  SLT target  = {v['SLT_target_f']:.3f} F")


if __name__ == "__main__":
    main()
