"""Figure 1: supercritical vs subcritical Hopf bifurcation (the Section 1.2 sketch).

POSTER VERSION -- sized for the CADEM0011 A1 e-poster.

The panel on the poster is 191 x 77 mm, so figsize is set to exactly that.
Every fontsize below is therefore a TRUE point size on the printed sheet.
Insert the PDF in PowerPoint, place at 100%, and never drag the corner
handles -- resizing after placement is what shrinks the text.

No data dependencies. Pure schematic.

Left  : supercritical -- stable cycles born at the Hopf point, nothing below.
Right : subcritical   -- stable and unstable branches meet at the fold with a
        vertical tangent; the unstable branch runs from the fold down to the
        Hopf point, spanning the bistable band. The envelope ends at the fold.

Writes: results/F1_schematic.pdf   (vector -- this is the one to place)
        results/F1_schematic.png   (300 dpi preview only)

Run:    python fig_schematic_poster.py
   or:  PYTHONPATH=src python examples/fig_schematic_poster.py
"""
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

# --------------------------------------------------------------------------
# output location
# --------------------------------------------------------------------------
OUTDIR = "results"
os.makedirs(OUTDIR, exist_ok=True)

# --------------------------------------------------------------------------
# poster geometry -- the panel this figure occupies on the A1 sheet
# --------------------------------------------------------------------------
MM = 25.4
PANEL_W_MM, PANEL_H_MM = 191.0, 77.0
FIGW, FIGH = PANEL_W_MM / MM, PANEL_H_MM / MM       # 7.52 x 3.03 in

# --------------------------------------------------------------------------
# typography -- Arial to match the poster, with a graceful fallback
# --------------------------------------------------------------------------
_available = {f.name for f in font_manager.fontManager.ttflist}
if "Arial" in _available:
    FONT = "Arial"
else:
    FONT = "DejaVu Sans"
    print("note: Arial not found, falling back to DejaVu Sans. "
          "Re-run on the machine that builds the poster for an exact match.")

plt.rcParams.update({
    "font.family":     FONT,
    "font.size":       15,
    "axes.titlesize":  17,
    "axes.labelsize":  17,
    "axes.linewidth":  1.2,
    "pdf.fonttype":    42,      # embed as TrueType, not Type 3
    "svg.fonttype":    "none",
})

BLUE, RED = "tab:blue", "tab:red"

# --------------------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(FIGW, FIGH), sharey=True)

# ---- (a) supercritical ---------------------------------------------------
a = ax[0]
Uh = 1.0
U = np.linspace(Uh, 1.9, 200)

a.plot([0.2, Uh], [0, 0], "k-", lw=2.6)                 # stable equilibrium
a.plot([Uh, 1.95], [0, 0], "k--", lw=1.8)               # unstable equilibrium
a.plot(U, 0.85 * np.sqrt(U - Uh), "-", color=BLUE, lw=3.0)
a.plot([Uh], [0], "ko", ms=10, mfc="none", mew=2.4)

a.annotate("flutter", (Uh, 0), textcoords="offset points",
           xytext=(0, -30), fontsize=15, ha="center")
a.annotate("stable LCO", (1.90, 0.85 * np.sqrt(0.90)), color=BLUE,
           textcoords="offset points", xytext=(0, 12), fontsize=15,
           ha="right")
a.set_title("Supercritical")

# ---- (b) subcritical -----------------------------------------------------
b = ax[1]
Uf, Uh2, Af = 0.55, 1.0, 0.42

b.axvspan(Uf, Uh2, color=RED, alpha=0.10)               # bistable band
b.plot([0.2, Uh2], [0, 0], "k-", lw=2.6)
b.plot([Uh2, 1.95], [0, 0], "k--", lw=1.8)

# unstable branch: fold -> Hopf, vertical tangent at the fold
Uu = np.linspace(Uf, Uh2, 200)
b.plot(Uu, Af * np.sqrt((Uh2 - Uu) / (Uh2 - Uf)), "--", color=RED, lw=2.4)

# stable branch: from the SAME fold point, vertical tangent
Us = np.linspace(Uf, 1.9, 200)
b.plot(Us, Af + 0.55 * np.sqrt(Us - Uf), "-", color=BLUE, lw=3.0)

b.plot([Uf], [Af], "v", color=BLUE, ms=14, zorder=5)
b.plot([Uh2], [0], "ko", ms=10, mfc="none", mew=2.4)

b.annotate("fold", (Uf, Af), textcoords="offset points",
           xytext=(-4, -34), fontsize=15, color=BLUE, ha="center")
b.annotate("flutter", (Uh2, 0), textcoords="offset points",
           xytext=(0, -30), fontsize=15, ha="center")
b.annotate("bistable", (0.5 * (Uf + Uh2), 0.86), fontsize=15,
           color=RED, ha="center")
b.annotate("unstable", (0.72, 0.30), fontsize=15, color=RED)
b.annotate("stable LCO", (1.90, Af + 0.55 * np.sqrt(1.35)), color=BLUE,
           textcoords="offset points", xytext=(0, 12), fontsize=15,
           ha="right")
b.set_title("Subcritical")

# ---- shared axis treatment ----------------------------------------------
for a_ in ax:
    a_.set_xlabel("airspeed")
    a_.set_xticks([])
    a_.set_yticks([])
    a_.set_xlim(0.2, 1.95)
    a_.set_ylim(-0.30, 1.25)
    a_.spines["top"].set_visible(False)
    a_.spines["right"].set_visible(False)

ax[0].set_ylabel("LCO amplitude")

plt.tight_layout(pad=0.4)

pdf = os.path.join(OUTDIR, "F1_schematic.pdf")
png = os.path.join(OUTDIR, "F1_schematic.png")
plt.savefig(pdf)
plt.savefig(png, dpi=300)
plt.close(fig)

print(f"-> {pdf}   ({PANEL_W_MM:.0f} x {PANEL_H_MM:.0f} mm, place at 100%)")
print(f"-> {png}   (preview only)")