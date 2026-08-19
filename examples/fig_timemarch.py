"""Figure: bistability in raw time histories (Michigan, ONERA rung, 10 m/s).

Inside the constructed bistable window (fold 9.594 < 10.0 < Hopf 10.875):
a 0.5 deg disturbance decays to the stable equilibrium while an 8 deg
disturbance is captured by the ~19 deg limit cycle. Two disturbances, two
fates -- the operational meaning of the fold, visible without any
bifurcation machinery.

Beware: the cubic must be set explicitly (beta = 1.326);
section_from_params() gives beta = 0.

Writes: results/F_timemarch.png
Run:    PYTHONPATH=src python examples/fig_timemarch.py   (~1-2 min)
"""
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

from mflco.model.params import TypicalSectionParameters
from mflco.model.michigan_params import calibrate_michigan, structural_zeta
from mflco.model.eom import structural_rhs
from mflco.aero.onera_stall import ONERAStall

cal = calibrate_michigan(zeta=structural_zeta())
p = TypicalSectionParameters(
    a=cal.a, x_alpha=cal.x_alpha, r_alpha_sq=cal.r_alpha_sq,
    omega_ratio=cal.omega_ratio, mu=cal.mu, beta=1.326, delta=0.0,
    zeta_h=cal.zeta, zeta_alpha=cal.zeta)
aero = ONERAStall(p, 0.0, N=6)
us = float(cal.ms_to_ustar(10.0))

fig, ax = plt.subplots(1, 2, figsize=(10.4, 3.5), sharex=True)
for axi, kick, ttl in (
        (ax[0], 0.5, "0.5 deg kick: decays -- equilibrium stable"),
        (ax[1], 8.0, "8 deg kick: captured by the ~19 deg LCO")):
    y0 = np.zeros(4 + aero.n_aero_states)
    y0[1] = np.radians(kick)
    sol = solve_ivp(structural_rhs, (0, 800), y0, args=(p, aero, us),
                    method="LSODA", rtol=1e-7, atol=1e-9,
                    t_eval=np.linspace(0, 800, 4000))
    axi.plot(sol.t, np.degrees(sol.y[1]), lw=0.8, color="tab:purple")
    axi.set_title(ttl, fontsize=9.5)
    axi.set_xlabel("non-dimensional time tau")
    axi.grid(alpha=0.25)
ax[0].set_ylabel("pitch [deg]")
fig.suptitle("Michigan wing at U = 10.0 m/s (ONERA rung): two disturbances, "
             "two fates -- bistability in the raw time histories",
             fontsize=10)
plt.tight_layout()
plt.savefig("results/F_timemarch.png", dpi=170)
print("-> results/F_timemarch.png")