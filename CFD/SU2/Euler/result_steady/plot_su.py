"""plot_su.py -- plot SU2 steady-Euler results (convergence + surface Cp).

Just hit the Run button (or `python plot_su.py`). Works from any folder, as
long as history.csv and surface_flow.csv sit in the SAME folder as this script.
Produces su2_plots.png next to the script.

Cp is computed from the conservative variables in surface_flow.csv
(Density, Momentum, Energy), since the steady run writes those, not pressure.
"""

import csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

HERE = Path(__file__).parent          # folder this script lives in -> finds the CSVs

# --- freestream (must match the .cfg) ---
GAMMA, R = 1.4, 287.058
P_INF, T_INF, MACH = 101325.0, 288.15, 0.3
RHO_INF = P_INF / (R * T_INF)
A_INF = np.sqrt(GAMMA * R * T_INF)
V_INF = MACH * A_INF
Q_INF = 0.5 * RHO_INF * V_INF ** 2          # dynamic pressure (= RefForce in history)

# --- history.csv ---
with open(HERE / "history.csv") as f:
    r = csv.reader(f)
    hdr = [h.strip().strip('"') for h in next(r)]
    H = np.array([[float(x) for x in row] for row in r if row])
col = {h: i for i, h in enumerate(hdr)}
it, rho = H[:, col["Inner_Iter"]], H[:, col["rms[Rho]"]]
CL, CD = H[:, col["CL"]], H[:, col["CD"]]
print(f"final  CL={CL[-1]:.4f}  CD={CD[-1]:.5f}  rms[Rho]={rho[-1]:.2f}")

# --- surface_flow.csv -> Cp ---
S = np.genfromtxt(HERE / "surface_flow.csv", delimiter=",", names=True)
x, y = S["x"], S["y"]
p = (GAMMA - 1) * (S["Energy"] - 0.5 * (S["Momentum_x"]**2 + S["Momentum_y"]**2) / S["Density"])
Cp = (p - P_INF) / Q_INF
up, lo = y >= 0, y < 0

# --- plot ---
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
ax[0].plot(it, rho, lw=1.2, color="steelblue")
ax[0].set(xlabel="iteration", ylabel="rms[Rho] (log10 residual)", title="Convergence history")
ax[0].grid(alpha=.3)
axc = ax[0].twinx(); axc.plot(it, CL, lw=1, color="darkorange", alpha=.8)
axc.set_ylabel("CL", color="darkorange")

ax[1].plot(x[up], Cp[up], ".", ms=3, color="crimson", label="upper")
ax[1].plot(x[lo], Cp[lo], ".", ms=3, color="steelblue", label="lower")
ax[1].invert_yaxis()
ax[1].set(xlabel="x/c", ylabel="Cp", title=f"Surface pressure (NACA0020, M={MACH}, 2 deg)")
ax[1].legend(); ax[1].grid(alpha=.3)

plt.tight_layout()
plt.savefig(HERE / "su2_plots.png", dpi=150)
print(f"saved {HERE / 'su2_plots.png'}")