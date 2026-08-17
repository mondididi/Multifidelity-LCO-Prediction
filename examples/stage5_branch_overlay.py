"""Stage 5c -- the two-wing overlay: branches, three-regime zero lines,
Hopf points, and the digitised experimental baselines.

Pure composition from the stage-5a/5b records plus examples/exp_digitised.py
(the single source of the experimental points). Zero lines carry the
subcritical picture for EVERY rung: solid to its fold (equilibrium the only
attractor), dotted fold->Hopf (bistable), dashed past the Hopf (equilibrium
unstable); supercritical rungs have no dotted segment, and Michigan's stall
rung is solid throughout (no Hopf at all). Rung zero-lines are stacked just
below zero so they stay readable; the experiment sits at true zero.

Bristol's experiment is HEAVE in mm (their diagram's DoF) on a twin right
axis: the velocities -- fold and Hopf -- are the like-for-like comparison.
Michigan's ONERA display is truncated at 26 deg (polar validity). The red
star is the beta amplitude anchor -- calibrated, not predicted.

Writes:  results/F_branch_overlay.png
Run:     PYTHONPATH=src python examples/stage5_branch_overlay.py   (seconds)
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from exp_digitised import draw_experiment

C = {"QS": "tab:blue", "QS_Stall": "tab:green", "Peters": "tab:orange",
     "ONERA": "tab:purple"}
BR = json.load(open("results/bristol_branches.json"))
MI = json.load(open("results/michigan_branches.json"))




fig, ax = plt.subplots(1, 2, figsize=(11.8, 4.6))

# ---- Bristol panel ----------------------------------------------------------------
b = ax[0]
for i, (key, lab) in enumerate((("QS", "QS"),
                                ("QS_Stall", "QS+Stall (schematic)"),
                                ("Peters", "Peters N=6"),
                                ("ONERA", "ONERA (schematic)"))):
    c = C[key]
    u, a = zip(*sorted(map(tuple, BR[key]["branch"])))
    b.plot(u, a, "o-", ms=3, lw=1.4, color=c, label=lab)
    b.plot([u[0]], [a[0]], "v", ms=8, color=c)
draw_experiment(b, "bristol", x_max=26.0, show_branch=False,
                label_prefix="experiment: eq / fold x / Hopf * "
                             "(amplitudes heave-mm, Fig. 6)")
b.set_title("Bristol wing: every rung folds (14.7-17.3); every Hopf within "
            "23.6-24.9", fontsize=9.5)
b.set_xlabel("airspeed U [m/s]")
b.set_ylabel("model pitch LCO amplitude [deg]")
b.legend(fontsize=7.3, loc="upper left")
b.grid(alpha=0.25)
b.set_xlim(13, 26)
b.set_ylim(-1.4, 31)

# ---- Michigan panel ---------------------------------------------------------------
m = ax[1]
for i, (key, lab) in enumerate((("QS", "QS (supercritical)"),
                                ("Peters", "Peters N=6 (supercritical)"),
                                ("ONERA",
                                 "ONERA, published laws (subcritical)"))):
    c = C[key]
    pts = sorted(map(tuple, MI[key]["branch"]))
    if key == "ONERA":
        pts = [t for t in pts if t[1] <= 26.0]
    u, a = zip(*pts)
    m.plot(u, a, "o-", ms=3, lw=1.4, color=c, label=lab)
    if key == "ONERA":
        m.plot([MI[key]["U_fold"]], [MI[key]["amp_at_fold"]], "v", ms=8,
               color=c)
        m.annotate("display truncated at polar validity (~25 deg)",
                   (u[-1], a[-1]), textcoords="offset points",
                   xytext=(-150, 8), fontsize=6.8, color=c)
draw_experiment(m, "michigan", x_max=14.6, label_prefix="experiment")
m.text(10.6, 28.0, "QS+Stall: equilibrium stable everywhere\n(no Hopf, no "
       "LCO -- tested null; green line)", fontsize=7.3, color=C["QS_Stall"])
m.set_title("Michigan wing: no attached rung folds; only dynamic stall is "
            "subcritical", fontsize=9.5)
m.set_xlabel("airspeed U [m/s]")
m.set_ylabel("pitch LCO amplitude [deg]")
m.legend(fontsize=7.3, loc="center left")
m.grid(alpha=0.25)
m.set_xlim(8.5, 14.6)
m.set_ylim(-1.4, 31)

fig.suptitle("Every executed rung on both wings: branches, three-regime "
             "equilibrium lines, Hopf points, digitised experiment",
             fontsize=11)
plt.tight_layout()
plt.savefig("results/F_branch_overlay.png", dpi=170)
print("-> results/F_branch_overlay.png")