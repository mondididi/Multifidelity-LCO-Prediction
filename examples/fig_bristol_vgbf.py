"""Figure: Bristol wing frequency/damping vs airspeed (QS and Peters).

The linearised story behind Section 4.1: the pitch and plunge branches
coalesce as airspeed rises -- a pattern set by the structural modes and
common to every rung -- and the damping crossing (flutter) sits within
+/-4% across the whole ladder on this wing, while the folds (drawn from
the committed stage-5 record) span 15%. Wake fidelity barely touches this
wing's flutter and moves its fold: the inverse of the Michigan wing.

Reads:  results/bristol_branches.json (fold/Hopf markers)
Writes: results/F_bristol_vgbf.png
Run:    PYTHONPATH=src python examples/fig_bristol_vgbf.py   (~1 min)
"""
import json

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from mflco.model.barton_params import BartonCal
from mflco.model.analysis import modal_analysis
from mflco.aero.quasi_steady import QuasiSteady
from mflco.aero.peters_finite import PetersFinite

cal = BartonCal(2)
p = cal.section()
REC = json.load(open("results/bristol_branches.json"))

fig, ax = plt.subplots(1, 2, figsize=(10.6, 4.1))
for name, aero, col in (("QS", QuasiSteady(p, 0.0), "tab:blue"),
                        ("Peters N=6", PetersFinite(p, 0.0, N=6),
                         "tab:orange")):
    f, d, ums = [], [], []
    for u in np.linspace(0.05, cal.ms_to_ustar(29.0), 420):
        try:
            fr, dp = modal_analysis(p, 0.0, u, aero)
        except ValueError:
            continue
        f.append(np.sort(fr)[:2])
        d.append(np.sort(dp)[:2])
        ums.append(cal.ustar_to_ms(u))
    f, d, ums = np.array(f), np.array(d), np.array(ums)
    hz = f * cal.omega_alpha / (2 * np.pi)
    for j in range(2):
        ax[0].plot(ums, hz[:, j], color=col, lw=1.4,
                   label=name if j == 0 else None)
        ax[1].plot(ums, d[:, j], color=col, lw=1.4,
                   label=name if j == 0 else None)
# for key, col in (("QS", "tab:blue"), ("Peters", "tab:orange")):
#     for a in ax:
#         a.axvline(REC[key]["hopf"], color=col, ls="--", lw=1.0, alpha=0.7)
#         a.axvline(REC[key]["U_fold"], color=col, ls=":", lw=1.2, alpha=0.7)
ax[1].axhline(0, color="0.4", lw=0.9)
for a in ax:
    a.set_xlabel("airspeed U [m/s]")
    a.grid(alpha=0.25)
    a.set_xlim(0, 29)
ax[0].set_ylabel("frequency [Hz]")
ax[1].set_ylabel("damping ratio")
ax[0].legend(fontsize=8)
ax[1].legend(fontsize=8)
fig.suptitle("Bristol Rig", fontsize=10.5)
plt.tight_layout()
plt.savefig("results/Fx_bristol_vgbf.png", dpi=170)
print("-> results/Fx_bristol_vgbf.png")