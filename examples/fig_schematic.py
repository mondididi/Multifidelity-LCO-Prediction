"""Figure: super- vs subcritical Hopf bifurcation (the Section 1.2 sketch).

Pure schematic -- no data dependencies. Left: supercritical (stable cycles
born at the Hopf point; nothing below flutter). Right: subcritical -- drawn
with the correct saddle-node topology: the stable large-amplitude branch
and the unstable branch MEET at the fold with a vertical tangent, and the
unstable branch (the basin boundary) runs continuously from the fold down
to the Hopf point, spanning the whole bistable band. The safe envelope
ends at the fold, not at flutter.

Writes: results/F1_schematic.png
Run:    PYTHONPATH=src python examples/fig_schematic.py   (seconds)
"""
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.7), sharey=True)

# ---- supercritical ---------------------------------------------------------------
a = ax[0]
Uh = 1.0
U = np.linspace(Uh, 1.9, 120)
a.plot([0.2, Uh], [0, 0], "k-", lw=2, label="stable equilibrium")
a.plot([Uh, 1.95], [0, 0], "k--", lw=1.4, label="unstable equilibrium")
a.plot(U, 0.85 * np.sqrt(U - Uh), "-", color="tab:blue", lw=2.2,
       label="stable LCO branch")
a.plot([Uh], [0], "ko", ms=7, mfc="none", mew=1.8)
a.annotate("Hopf (flutter)", (Uh, 0), textcoords="offset points",
           xytext=(-10, -18), fontsize=9)
a.set_title("Supercritical: nothing below flutter", fontsize=10)
a.legend(fontsize=7, loc="upper left")

# ---- subcritical: correct saddle-node topology -----------------------------------
b = ax[1]
Uf, Uh2, Af = 0.55, 1.0, 0.42
b.axvspan(Uf, Uh2, color="tab:red", alpha=0.10)
b.plot([0.2, Uh2], [0, 0], "k-", lw=2)
b.plot([Uh2, 1.95], [0, 0], "k--", lw=1.4)
# unstable branch: fold -> Hopf, vertical tangent at the fold
Uu = np.linspace(Uf, Uh2, 160)
b.plot(Uu, Af * np.sqrt((Uh2 - Uu) / (Uh2 - Uf)), "--", color="tab:red",
       lw=1.9, label="unstable LCO (basin boundary)")
# stable branch: continues from the SAME fold point, vertical tangent
Us = np.linspace(Uf, 1.9, 160)
b.plot(Us, Af + 0.55 * np.sqrt(Us - Uf), "-", color="tab:blue", lw=2.2,
       label="stable LCO branch")
b.plot([Uf], [Af], "v", color="tab:blue", ms=11, zorder=5)
b.plot([Uh2], [0], "ko", ms=7, mfc="none", mew=1.8)
b.annotate("fold", (Uf, Af), textcoords="offset points", xytext=(-32, 2),
           fontsize=9)
b.annotate("Hopf (flutter)", (Uh2, 0), textcoords="offset points",
           xytext=(-10, -18), fontsize=9)
b.annotate("bistable", (0.5 * (Uf + Uh2), 0.68), fontsize=9,
           color="tab:red", ha="center")
b.set_title("Subcritical: the envelope ends at the fold", fontsize=10)
b.legend(fontsize=7, loc="upper left")

for a_ in ax:
    a_.set_xlabel("airspeed")
    a_.set_xticks([])
    a_.set_yticks([])
    a_.set_ylim(-0.12, 1.02)
ax[0].set_ylabel("LCO amplitude")
plt.tight_layout()
plt.savefig("results/F1_schematic.png", dpi=170)
print("-> results/F1_schematic.png")