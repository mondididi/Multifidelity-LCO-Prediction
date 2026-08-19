"""Figure: Michigan wing frequency/damping vs airspeed (QS and Peters).

The linearised story behind Section 5.1: the same coalescence pattern as
the Bristol wing -- set by the structural modes -- but here wake fidelity
moves the damping crossing enormously: QS flutters at 9.665 (-26.7%) while
Peters' six inflow states restore the lag and land at 13.148 (-0.32%,
unfitted) against the measured 13.19. Fold/Hopf markers are read from the
committed stage-5 record; the measured wind-off modes (5.30 / 6.20 Hz) are
overlaid at the left edge.

Reads:  results/michigan_branches.json
Writes: results/F_michigan_vgbf.png
Run:    PYTHONPATH=src python examples/fig_michigan_vgbf.py   (~1 min)
"""
import json

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from mflco.model.params import TypicalSectionParameters
from mflco.model.michigan_params import calibrate_michigan, structural_zeta
from mflco.model.analysis import modal_analysis
from mflco.aero.quasi_steady import QuasiSteady
from mflco.aero.peters_finite import PetersFinite

cal = calibrate_michigan(zeta=structural_zeta())
p = TypicalSectionParameters(
    a=cal.a, x_alpha=cal.x_alpha, r_alpha_sq=cal.r_alpha_sq,
    omega_ratio=cal.omega_ratio, mu=cal.mu, beta=1.326, delta=0.0,
    zeta_h=cal.zeta, zeta_alpha=cal.zeta)
REC = json.load(open("results/michigan_branches.json"))
MEASURED_HZ = (5.30, 6.20)           # wind-off modes (Garcia Perez)
EXP_FLUTTER = 13.19

fig, ax = plt.subplots(1, 2, figsize=(10.6, 4.1))
for name, aero, col in (("QS", QuasiSteady(p, 0.0), "tab:blue"),
                        ("Peters N=6", PetersFinite(p, 0.0, N=6),
                         "tab:orange")):
    f, d, ums = [], [], []
    for u in np.linspace(0.05, cal.ms_to_ustar(15.0), 420):
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
# for hzm in MEASURED_HZ:
#     ax[0].plot([0.15], [hzm], "k<", ms=7)
# ax[0].plot([], [], "k<", ms=7, label="measured wind-off modes")
for a in ax:
    #a.axvline(EXP_FLUTTER, color="r", ls=":", lw=1.4)
    a.set_xlabel("airspeed U [m/s]")
    a.grid(alpha=0.25)
    a.set_xlim(0, 15)
ax[1].axhline(0, color="0.4", lw=0.9)
ax[0].set_ylabel("frequency [Hz]")
ax[1].set_ylabel("damping ratio")
ax[0].legend(fontsize=8)
ax[1].legend(fontsize=8)
fig.suptitle("Michigan rig", fontsize=10)
plt.tight_layout()
plt.savefig("results/Fx_michigan_vgbf.png", dpi=170)
print("-> results/Fx_michigan_vgbf.png")