# R-22 Superheat Charging Assistant
 
A single-file web tool for charging fixed-orifice R-22 systems using the superheat method. Log gauge readings as you go, and it tells you whether to add refrigerant or stop, tracks the trend across a session, and remembers history per site so estimates improve over time.
 
**Live tool:** https://omouraenko.github.io/r22-tool/
 
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
