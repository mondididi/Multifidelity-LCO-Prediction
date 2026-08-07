"""Stage 2b -- Michigan delta sweep by fold-velocity down-sweep (the crux).

Supersedes the contique route of stage2_delta_sweep.py for the large-delta
cases: continuation hangs when Newton's trial states hit the violently stiff
cubic at large amplitude; the time-marching down-sweep has no Newton and
cannot hang. Cross-checks: delta = 0 and 3.33 reproduce the continuation
verdicts (supercritical); the detector itself is verified at stage 0 against
Barton's published model (Hopf to 0.08%, fold to 0.6%).

Protocol per delta: beta re-calibrated to the amplitude anchor (14 deg at
13.5 m/s) so the fold stays a PREDICTION; b1 fixed so the Hopf stays
13.148 m/s at every delta; then the down-sweep reads the fold directly.
Anchor amplitude measure: HALF PEAK-TO-PEAK of the late window -- robust to
the asymmetric orbits large delta produces (the mean-positive-peak measure of
stage2 becomes unreachable there; for near-symmetric small-delta orbits the
two measures coincide and verdicts are identical). delta = 0 keeps the
stage-2 calibrated beta = 1.326 for continuity.

Result (first run 7/8/2026): supercritical at delta = 0, 3.33, 6.67; at the
physical bracket edge delta = 10 a genuine but tiny fold appears -- window
0.02-0.05 m/s (long-march verified ~857 periods at 13.125; +/-20 deg kicks at
12.5 and 11.9 m/s all decay) vs the rig's 1.34 m/s. Claim to carry: no
attached-flow model with any admissible member of the rig's documented
cubic-stiffness family reproduces more than ~4% of the observed bistable
window. Report the delta = 10 window as a BOUND at the 0.05 m/s resolution.

Cold run ~4 min (beta calibrations dominate); cached incrementally to the
system temp dir -- rerun plots instantly, delete the JSON to recompute.

"""
import json
import os
import tempfile

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from mflco.model.params import TypicalSectionParameters
from mflco.model.michigan_params import calibrate_michigan, structural_zeta
from mflco.model.eom import structural_rhs
from mflco.model.fold_detector import classify_orbit, down_sweep
from mflco.aero.peters_finite import PetersFinite

U_HOPF_MS   = 13.148   # this build's linear flutter; delta-invariant (b1 fixed)
RIG_FLUTTER = 13.19    # Garcia Perez J063736
RIG_FOLD    = 11.85    # Garcia Perez turning point
N_INFLOW    = 6        # Peters states, as stages 1-2
U_ANCHOR    = 13.5     # amplitude anchor speed, above flutter
AMP_TARGET  = 14.0     # deg, half peak-to-peak anchor
BETA_D0     = 1.326    # stage-2 calibrated beta at delta = 0 (kept: measures
                       # coincide for the symmetric orbit)
DELTAS_DEG  = [0.0, 3.33, 6.67, 10.0]
CACHE = os.path.join(tempfile.gettempdir(), "stage2b_delta_downsweep.json")
PNG   = os.path.join(tempfile.gettempdir(), "stage2b_delta_downsweep.png")

cal = calibrate_michigan(zeta=structural_zeta())


def _section(beta, delta_rad):
    return TypicalSectionParameters(
        a=cal.a, x_alpha=cal.x_alpha, r_alpha_sq=cal.r_alpha_sq,
        omega_ratio=cal.omega_ratio, mu=cal.mu, beta=beta, delta=delta_rad,
        zeta_h=cal.zeta, zeta_alpha=cal.zeta)


def _blowup(tau, y, *a):                     # calibration-only guard
    return abs(y[1]) - np.radians(80.0)
_blowup.terminal = True


def _settle_amp(beta, delta_rad, U_ms, tau_end=650.0):
    """Half peak-to-peak pitch [deg] of the late window at U_ms."""
    p = _section(beta, delta_rad)
    aero = PetersFinite(p, 0.0, N=N_INFLOW)
    y0 = np.zeros(4 + N_INFLOW)
    y0[1] = np.radians(10.0)
    sol = solve_ivp(structural_rhs, (0.0, tau_end), y0,
                    args=(p, aero, float(cal.ms_to_ustar(U_ms))),
                    method="LSODA", rtol=1e-7, atol=1e-9, events=_blowup,
                    t_eval=np.linspace(0.0, tau_end, 16000))
    if sol.status == 1:
        return 90.0                          # blew past 80 deg: "amplitude huge"
    late = sol.y[1][sol.t > 0.85 * tau_end]
    return float(np.degrees(0.5 * (late.max() - late.min())))


def _load():
    try:
        return json.load(open(CACHE))
    except (FileNotFoundError, json.JSONDecodeError):
        return {}


def _save(c):
    json.dump(c, open(CACHE, "w"), indent=1)


# --- compute (load-or-run, incremental) ------------------------------------------
cache = _load()
for d_deg in DELTAS_DEG:
    key = f"{d_deg:.2f}"
    c = cache.setdefault(key, {})
    if "beta" not in c:
        if d_deg == 0.0:
            c["beta"] = BETA_D0
        else:
            print(f"[delta={d_deg}] calibrating beta to {AMP_TARGET} deg @ "
                  f"{U_ANCHOR} m/s ...", flush=True)
            c["beta"] = round(brentq(
                lambda b: _settle_amp(b, np.radians(d_deg), U_ANCHOR)
                - AMP_TARGET, 0.05, 12.0, xtol=0.02), 3)
        print(f"[delta={d_deg}] beta = {c['beta']}", flush=True)
        _save(cache)
    if "branch" not in c:
        print(f"[delta={d_deg}] down-sweep ...", flush=True)
        p = _section(c["beta"], np.radians(d_deg))
        aero = PetersFinite(p, 0.0, N=N_INFLOW)
        y0 = np.zeros(4 + N_INFLOW)
        y0[1] = np.radians(10.0)
        argU = lambda U: (p, aero, float(cal.ms_to_ustar(U)))
        br, U_fold, amp_f, secs = down_sweep(
            structural_rhs, argU, y0, U_ANCHOR, 0.1, 10.5, 0.05)
        window = U_HOPF_MS - U_fold
        c.update(branch=[(round(u, 3), round(a, 2)) for u, a in br],
                 U_fold=round(U_fold, 3), amp_at_fold=round(amp_f, 2),
                 window=round(window, 3), seconds=round(secs, 1),
                 verdict=("SUBCRITICAL" if (window > 0.05 and amp_f > 1.5)
                          else "SUPERCRITICAL"))
        print(f"[delta={d_deg}] U_fold {c['U_fold']}, window {c['window']}, "
              f"{c['verdict']} [{c['seconds']}s]", flush=True)
        _save(cache)

# --- report + figure --------------------------------------------------------------
print(f"\n{'delta':>6} {'beta':>7} {'U_fold':>7} {'amp@fold':>8} {'window':>7} "
      f"{'cost[s]':>8}  verdict")
for d_deg in DELTAS_DEG:
    c = cache[f"{d_deg:.2f}"]
    print(f"{d_deg:>6.2f} {c['beta']:>7.3f} {c['U_fold']:>7.3f} "
          f"{c['amp_at_fold']:>8.2f} {c['window']:>7.3f} {c['seconds']:>8.1f}"
          f"  {c['verdict']}")
print(f"   rig {'--':>7} {RIG_FOLD:>7.2f} {'~4':>8} {1.34:>7.2f} {'--':>8}"
      f"  subcritical")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(8.4, 5.2))
for d_deg in DELTAS_DEG:
    c = cache[f"{d_deg:.2f}"]
    br = np.array(sorted(c["branch"]))
    ax.plot(br[:, 0], br[:, 1], "o-", ms=3,
            label=f"delta={d_deg:.2f}  beta={c['beta']:.2f}  "
                  f"({c['verdict'][:5]})")
ax.axvline(RIG_FLUTTER, color="r", ls="--", lw=1.1, label="rig flutter 13.19")
ax.axvline(RIG_FOLD, color="r", ls=":", lw=1.4, label="rig fold 11.85")
ax.axvline(U_HOPF_MS, color="0.5", ls="-.", lw=1.0, label="model Hopf 13.148")
ax.set_xlabel("airspeed U [m/s]")
ax.set_ylabel("pitch LCO amplitude [deg]")
ax.set_title("Stage 2b: does the trim quadratic produce the fold? "
             "(down-sweep, all delta)")
ax.legend(fontsize=8)
ax.grid(alpha=0.25)
plt.tight_layout()
plt.savefig(PNG, dpi=170)
print(f"\ncache  -> {CACHE}")
print(f"figure -> {PNG}")