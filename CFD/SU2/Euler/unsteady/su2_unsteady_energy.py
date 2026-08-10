"""Post-process the SU2 prescribed-pitch campaign -> F8, the inviscid ceiling.

Reads the unsteady history CSVs (one row per physical time step), reconstructs
the prescribed incidence alpha(t) = AOA + AMPL*sin(omega*t), and builds the
Cl-alpha loops from the LAST TWO periods (the first two absorb the start-up
transient). Two overlays make the argument:

  small (2 deg):  the mflco Peters N=6 model driven by the SAME prescribed
                  motion -- if the Euler loop lands on the Peters loop, the
                  inviscid ceiling and the finite-state model are the same
                  physics at the rig's reduced frequency, executed two ways.
  large (14 deg): the STATIC polar band swept by the motion (-4..+24 deg) --
                  the Euler loop stays an attached-flow ellipse straight
                  through the range where the real airfoil is deep in stall:
                  no stall hysteresis, no fold mechanism, at ANY inviscid
                  fidelity. That is the ceiling statement as a picture.

Also prints the cycle work coefficient W = closed-integral Cm dalpha (about
the pitch axis) per cycle: negative = aerodynamic damping (energy extracted
from the motion), the attached-flow expectation at k = 0.286.

Peters overlay: the structural states are PRESCRIBED (pure pitch about the
elastic axis, a = -0.5), only the six inflow states are integrated, and Cl is
recovered from the generalized force via Cl = -Q_xi * pi * mu / U*^2 -- the
exact inverse of the lift scaling in the force routines.

Run (from the folder holding history_small.csv / history_large.csv):
  PYTHONPATH=<repo>\\src python su2_unsteady_energy.py [case_dir]
Outputs F8_euler_loops.png + su2_unsteady_energy.json next to the histories.
"""
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- campaign constants (must match the cfg headers) -------------------------------
AOA_DEG    = 10.0          # mean incidence: the rig trim
OMEGA      = 58.42         # rad/s in CFD units -> k = 0.286 at M = 0.3, c = 1 m
CASES      = {"small": 2.0, "large": 14.0}     # pitching amplitudes [deg]
N_LAST     = 2             # periods used for the loops (of 4 simulated)
A_EA       = -0.5          # elastic axis in semichords: pitch axis = c/4

CASE_DIR = sys.argv[1] if len(sys.argv) > 1 else "."


def _read_history(path):
    """SU2 history CSV: quoted, space-padded headers; return (t, CL, CMz)."""
    import csv as _csv
    with open(path, newline="") as f:
        rows = list(_csv.reader(f))
    hdr = [h.strip().strip('"').strip() for h in rows[0]]
    data = np.array([[float(x) for x in r] for r in rows[1:] if r])
    col = {name: hdr.index(name) for name in hdr}
    t = data[:, col["Cur_Time"]] if "Cur_Time" in col else None
    cl = data[:, col["CL"]]
    cm = data[:, col["CMz"]] if "CMz" in col else data[:, col["CM"]]
    return t, cl, cm


def _last_periods(t, alpha, y, n=N_LAST):
    T = 2 * np.pi / OMEGA
    m = t >= t[-1] - n * T
    return alpha[m], y[m], m


def _cycle_work(alpha_rad, cm):
    """W = closed-integral Cm d(alpha) over the retained window / N_LAST."""
    return float(np.trapz(cm, alpha_rad) / N_LAST)


# --- Peters overlay: same motion, inflow states only -------------------------------
def peters_prescribed_cl(amp_deg, n_per=2, pts=600):
    from scipy.integrate import solve_ivp
    from mflco.model.michigan_params import (calibrate_michigan,
                                             section_from_params,
                                             structural_zeta)
    from mflco.aero.peters_finite import PetersFinite
    cal = calibrate_michigan(zeta=structural_zeta())
    p = section_from_params(cal)
    aero = PetersFinite(p, 0.0, N=6)
    # prescribed pitch about the EA in package time tau = omega_alpha * t:
    # match the REDUCED frequency: omega_tau chosen so k = 0.286 at U*_anchor
    U_star = float(cal.ms_to_ustar(12.5))
    k = 0.286
    om_tau = k * U_star                       # alpha = A sin(om_tau * tau)
    A = np.radians(amp_deg)

    def ys(tau):
        al = A * np.sin(om_tau * tau)
        ald = A * om_tau * np.cos(om_tau * tau)
        return np.array([0.0, al, 0.0, ald])

    def rhs(tau, lam):
        return aero.aero_rhs(tau, ys(tau), lam, U_star)

    T = 2 * np.pi / om_tau
    sol = solve_ivp(rhs, (0.0, (n_per + 2) * T), np.zeros(6), method="LSODA",
                    rtol=1e-8, atol=1e-10,
                    t_eval=np.linspace(2 * T, (n_per + 2) * T, pts))
    cl = []
    for tau, lam in zip(sol.t, sol.y.T):
        Q = aero.forces(tau, ys(tau), lam, U_star)
        cl.append(-Q[0] * np.pi * p.mu / U_star ** 2)   # invert lift scaling
    al = np.degrees(A * np.sin(om_tau * sol.t))
    return al, np.asarray(cl)


# --- static polar band for the large case ------------------------------------------
def polar_band():
    from mflco.aero.qs_stall import load_polar
    ag, cl = load_polar()
    m = (ag >= -8) & (ag <= 28)
    return ag[m], cl[m]


# --- main --------------------------------------------------------------------------
out = {}
fig, ax = plt.subplots(1, 2, figsize=(11.4, 4.8))
for j, (case, amp) in enumerate(CASES.items()):
    path = os.path.join(CASE_DIR, f"history_{case}.csv")
    if not os.path.exists(path):
        ax[j].set_title(f"{case}: history_{case}.csv not found")
        continue
    t, cl, cm = _read_history(path)
    if t is None:                             # fall back: uniform steps
        t = np.arange(len(cl)) * 8.4023e-4
    al = AOA_DEG + amp * np.sin(OMEGA * t)
    al_w, cl_w, m = _last_periods(t, al, cl)
    _, cm_w, _ = _last_periods(t, al, cm)
    W = _cycle_work(np.radians(al_w), cm_w)
    out[case] = dict(cycle_work_Cm_dalpha=W, n_steps=int(m.sum()))
    ax[j].plot(al_w, cl_w, "-", lw=1.6, color="tab:red",
               label=f"SU2 Euler, k=0.286 (last {N_LAST} periods)")
    ax[j].set_xlabel("incidence [deg]")
    ax[j].set_ylabel("Cl")
    ax[j].grid(alpha=0.25)
    ax[j].set_title(f"{case}: {amp:.0f} deg about trim   "
                    f"W_cyc = {W:+.4f} (Cm dalpha)")

# overlays
try:
    al_p, cl_p = peters_prescribed_cl(CASES["small"])
    # Peters Cl is a PERTURBATION about trim: shift onto absolute incidence
    ax[0].plot(AOA_DEG + al_p, cl_p + float(np.interp(0, al_p, cl_p) * 0 + 0),
               "--", lw=1.4, color="tab:orange",
               label="Peters N=6, same motion (perturbation, shifted to trim)")
    ax[0].legend(fontsize=8)
except Exception as e:                        # mflco not on path: still plot SU2
    ax[0].legend(fontsize=8)
    print(f"[peters overlay skipped: {e}]")
try:
    ag, cls = polar_band()
    ax[1].plot(ag, cls, ":", lw=1.6, color="tab:green",
               label="static polar (what the real airfoil does)")
    ax[1].legend(fontsize=8)
except Exception as e:
    ax[1].legend(fontsize=8)
    print(f"[polar overlay skipped: {e}]")

fig.suptitle("F8  The inviscid ceiling: attached-flow loops straight through "
             "the stall range", fontsize=11)
plt.tight_layout()
png = os.path.join(CASE_DIR, "F8_euler_loops.png")
plt.savefig(png, dpi=170)
json.dump(out, open(os.path.join(CASE_DIR, "su2_unsteady_energy.json"), "w"),
          indent=1)
print(f"figure -> {png}")
print(json.dumps(out, indent=1))
