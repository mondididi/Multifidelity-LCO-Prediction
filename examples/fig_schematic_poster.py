"""Figure: super- vs subcritical Hopf bifurcation (the Section 1.2 sketch).

POSTER VERSION. The panel on the A1 sheet is 191 x 77 mm, so figsize is set
to exactly that. Every fontsize below is therefore a TRUE point size on the
printed poster. Place at 100% in PowerPoint and never drag the handles.

Legends replaced by direct labels on the curves: on a poster the eye should
not have to travel to a key and back.

Writes: results/F1_schematic.pdf  (vector - line art, so no dpi needed)
Run:    PYTHONPATH=src python examples/fig_schematic_poster.py
"""
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

MM = 25.4
FIGW, FIGH = 191 / MM, 77 / MM          # 7.52 x 3.03 in -> the placed size

plt.rcParams.update({
    "font.family": "Arial",             # matches the poster
    "font.size": 15,
    "axes.titlesize": 17,
    "axes.labelsize": 17,
    "axes.linewidth": 1.2,
})

BLUE, RED = "tab:blue", "tab:red"

fig, ax = plt.subplots(1, 2, figsize=(FIGW, FIGH), sharey=True)

# ---- supercritical ---------------------------------------------------------
a = ax[0]
Uh = 1.0
U = np.linspace(Uh, 1.9, 200)
a.plot([0.2, Uh], [0, 0], "k-", lw=2.6)
a.plot([Uh, 1.95], [0, 0], "k--", lw=1.8)
a.plot(U, 0.85 * np.sqrt(U - Uh), "-", color=BLUE, lw=3.0)
a.plot([Uh], [0], "ko", ms=10, mfc="none", mew=2.4)

a.annotate("flutter", (Uh, 0), textcoords="offset points", xytext=(0, -30),
           fontsize=15, ha="center")
a.annotate("stable LCO", (1.62, 0.85 * np.sqrt(0.62)), color=BLUE,
           textcoords="offset points", xytext=(-6, 12), fontsize=15)
a.set_title("Supercritical")

# ---- subcritical -----------------------------------------------------------
b = ax[1]
Uf, Uh2, Af = 0.55, 1.0, 0.42
b.axvspan(Uf, Uh2, color=RED, alpha=0.10)
b.plot([0.2, Uh2], [0, 0], "k-", lw=2.6)
b.plot([Uh2, 1.95], [0, 0], "k--", lw=1.8)

Uu = np.linspace(Uf, Uh2, 200)
b.plot(Uu, Af * np.sqrt((Uh2 - Uu) / (Uh2 - Uf)), "--", color=RED, lw=2.4)

Us = np.linspace(Uf, 1.9, 200)
b.plot(Us, Af + 0.55 * np.sqrt(Us - Uf), "-", color=BLUE, lw=3.0)

b.plot([Uf], [Af], "v", color=BLUE, ms=14, zorder=5)
b.plot([Uh2], [0], "ko", ms=10, mfc="none", mew=2.4)

b.annotate("fold", (Uf, Af), textcoords="offset points", xytext=(-4, -34),
           fontsize=15, color=BLUE, ha="center")
b.annotate("flutter", (Uh2, 0), textcoords="offset points", xytext=(0, -30),
           fontsize=15, ha="center")
b.annotate("bistable", (0.5 * (Uf + Uh2), 0.86), fontsize=15, color=RED,
           ha="center")
b.annotate("unstable", (0.72, 0.30), fontsize=15, color=RED)
b.annotate("stable LCO", (1.55, Af + 0.55 * np.sqrt(1.0)), color=BLUE,
           textcoords="offset points", xytext=(-14, 6), fontsize=15)
b.set_title("Subcritical")

for a_ in ax:
    a_.set_xlabel("airspeed")
    a_.set_xticks([])
    a_.set_yticks([])
    a_.set_ylim(-0.30, 1.25)
    a_.set_xlim(0.2, 1.95)
    for side in ("top", "right"):
        a_.spines[side].set_visible(False)

ax[0].set_ylabel("LCO amplitude")
plt.tight_layout(pad=0.4)
plt.savefig("results/F1_schematic.pdf")
plt.savefig("results/F1_schematic.png", dpi=300)   # preview only
print("-> results/F1_schematic.pdf")