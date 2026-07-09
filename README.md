# R-22 Superheat Charging Assistant

A single-file web tool for charging fixed-orifice R-22 systems using the superheat method. Log gauge readings as you go, and it tells you whether to add refrigerant or stop, tracks the trend across a session, and remembers history per site to estimates freon loss over time.

**Live tool:** <https://omouraenko.github.io/r22-tool/>

## What it does

- Walks through baseline → charging steps → results for a recharge, or a quick read-only checkup
- Computes target superheat, actual superheat, subcooling, and condensing split from your pressure/temperature readings
- Tells you to add refrigerant, stop, or flags a problem worth investigating
- Tracks refrigerant added per session (via jug weight) and estimates how much more is needed
- Keeps history per site (system), and uses that history to improve its estimates over time
- Charts session and multi-visit trends, including a rough "next recharge due" forecast based on estimated leak rate
- Everything is editable after the fact — you can correct or delete any past reading if you find an error
- Works offline once loaded; all data stays on your device (see below)

## Your data

This tool stores everything locally in your browser (`localStorage`) — nothing is sent to a server. That means:

- Your site/session history is tied to *this specific browser on this specific device*. It won't show up if you open the same link on a different phone or in a different browser.
- Clearing your browser's site data will erase it.
- Use the **Backup** button regularly (Site section, and the Backup panel at the bottom of the page) to export a `.json` file as a safety copy, and **Restore** to bring it back.

## Using tool on iPhone

The tool was tested on an iPhone. A few additional steps can improve its usability:

- Enable the Keep Screen Awake toggle at the top of the page to prevent the screen from going to sleep while the tool is in use.
- Create a Home Screen shortcut for quick access.
- Use the iPhone Shortcuts app to save clipboard data to a backup file.

### Add to Home Screen

The tool works like a regular app if you add it to your Home Screen — it opens full-screen without Safari's address bar and gets its own icon.

1. Open the [live tool](https://omouraenko.github.io/r22-tool/) in **Safari** (this only works in Safari, not Chrome or other browsers on iOS).
2. Tap the **Share** button (the square with an arrow pointing up) in the toolbar.
3. Scroll down and tap **Add to Home Screen**.
4. Tap **Add** in the top-right corner.

The icon will appear on your Home Screen like any other app. Note that your saved history is tied to Safari's storage for this site — if you clear Safari's site data, it'll erase the Home Screen app's data too, since they share the same storage.

### Saving backups into files using Shortcuts

iOS Safari can't save a file directly to disk, so when you tap **Backup**, the tool falls back to showing the JSON in a dialog with a **Copy to clipboard** button — from there you'd normally paste it into Notes or the Mail app. A **Shortcuts** automation turns that into a proper `.json` file in iCloud Drive with one tap, and avoids overwriting older backups by mistake.

#### Option A: install the ready-made shortcut

1. On your iPhone, open this link in **Safari**: [R22 backup.shortcut](https://raw.githubusercontent.com/omouraenko/r22-tool/main/R22%20backup.shortcut).
2. Safari downloads it — tap the **downloads** icon (circle with a down arrow) in the toolbar, then tap the downloaded file to hand it off to the Shortcuts app.
3. If iOS shows an **"Untrusted Shortcut"** warning (common for shortcuts downloaded outside of Apple's iCloud share links), go to **Settings → Shortcuts → Advanced → Allow Untrusted Shortcuts**, turn it on, then reopen the downloaded file.
4. Shortcuts shows you the list of actions — review them, then tap **Add Shortcut**.
5. Open the shortcut once in edit mode and re-pick your own destination folder in the **Get Contents of** and **Save File** actions — the original iCloud Drive folder reference is tied to the device it was built on and won't carry over to yours.

#### Option B: build it yourself

Build a new shortcut (Shortcuts app → **+**) named e.g. "R22 backup" with these actions, in order:

1. **Current Date**
2. **Format Date** — format the date from step 1 (e.g. `yyyy-MM-dd`).
3. **Text** — `backup-` followed by the **Formatted Date** variable.
4. **Set Variable** — name it **Name**, value = the Text from step 3 (e.g. `backup-2026-07-08`).
5. **Text** — the **Name** variable followed by `.json`.
6. **Set Variable** — name it **Full Name**, value = the Text from step 5 (e.g. `backup-2026-07-08.json`).
7. **Get Contents of** — pick the iCloud Drive folder you want backups saved to.
8. **Filter Files** — filter **Folder Contents** (the output of step 7) where *all* of: **Name is Name** and **File Extension is json**; turn on **Limit**, set to **Get 1 File**.
9. **If** — condition: **Files** (the Filter output) **has any value**:
   - **Choose from Menu** — prompt: "A file named Full Name already exists. Overwrite?" with options **Overwrite** and **Cancel**.
     - Under **Overwrite**: add **Delete Files**, input = **Files**.
     - Under **Cancel**: add **Stop This Shortcut**.
10. **Get Clipboard**
11. **Save File** — input = **Clipboard**, destination = the same folder as step 7.
12. **Rename File** — input = **Saved File** (from step 11), new name = **Full Name**.

To use it: in the tool, tap **Backup** → **Copy to clipboard** in the dialog, then switch over and run the "R22 backup" shortcut (from the Shortcuts app, a Home Screen icon, or the Share Sheet).

A few things that trip this up:

- **Save File** always saves as `.txt` when the input is plain text — that's why step 12 renames it afterward to force the `.json` extension.
- **Deleting an existing file still asks for one-time confirmation** even with the Delete Files action's own "Immediately Delete" option on — iOS gates it separately under **Settings → Shortcuts → Advanced → Allow Deleting Without Confirmation** (or similar wording, varies by iOS version).
- The destination folder picked in steps 7 and 11 is tied to your iCloud account/device — anyone installing the shortcut (via Option A or B) needs to point those two actions at their own folder before running it.

## Use at your own risk

This is a personal tool, not a certified instrument or substitute for proper training:

- **Charging refrigerant involves real hazards** — pressurized cylinders, frostbite from liquid refrigerant, and system damage from overcharging. Follow standard safety practice: charge vapor only, cylinder upright, compressor running, and stop if anything about the system or your readings seems off.
- **The math here is a planning aid, not a guarantee.** Target superheat, leak-rate forecasts, and "how much to add" estimates are all approximations based on the readings you enter — bad input (a slipped gauge, a stuck thermometer) gives a bad output, and the tool has no way to catch that.
- **Handling R-22 is regulated.** In the US, purchasing or working with R-22 requires EPA Section 608 certification. This tool doesn't check for or substitute for that — it assumes the person using it is already qualified to do the work.
- No warranty of any kind. Verify readings against your own judgment and training before acting on anything this tool tells you.

## Running it yourself

It's one self-contained HTML file — no build step, no dependencies, no backend. To update the live version, replace `index.html` in this repo with the new file and commit.

## License

MIT — see [LICENSE](LICENSE). Use, copy, and modify freely; it comes with no warranty (see the disclaimer above).
