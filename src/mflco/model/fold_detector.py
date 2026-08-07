"""Fold-velocity detector: time-marching down-sweep for LCO branches.

Reads the fold velocity directly: settle an LCO above the Hopf point, step the
airspeed DOWN, warm-starting each step from the previous orbit, and record the
speed at which the LCO collapses to equilibrium. That collapse speed is the
fold velocity -- the same protocol used experimentally (Garcia Perez Fig. 9
down-sweep; Tartaruga Fig. 6 decreasing-speed points), so numerical and
experimental numbers are like-for-like. No continuation, no Newton iteration:
the large-amplitude stiffness hang that stops arclength continuation cannot
occur here.

Classification per airspeed (chunked march, LSODA):
    LCO      -- chunk-over-chunk peak pitch agrees within rel_tol
    collapse -- peak pitch below floor_deg; or small (< 3 deg) and
                monotonically decaying (early exit); or still monotonically
                decaying after max_chunks (critical-slowing guard)

Verification provenance (stage 0): fed Barton's published model (AeroelasticCBC
src/model.jl case 2) through this detector via the mflco pipeline, recovers
Hopf 23.98 m/s and fold 16.92 m/s against the independent full-Sears port's
24.0 / 16.81 and the Tartaruga rig's experimental ~24 / ~16 m/s. Michigan
delta = 0 reproduces the stage-2 continuation verdict (supercritical).

State convention: pitch is state index 1, consistent with eom.structural_rhs
([xi, alpha, xi_dot, alpha_dot, aero...]). Resolution: fold bracketed to
`refine` m/s; sub/supercritical calls at the Hopf are meaningful only to that
resolution.
"""

import time

import numpy as np
from scipy.integrate import solve_ivp


# --- tuning defaults (reasons inline) --------------------------------------------
FLOOR_DEG  = 0.5     # below this peak pitch the orbit has collapsed
REL_TOL    = 0.02    # chunk-over-chunk amplitude agreement -> converged LCO
MAX_CHUNKS = 8       # then decide by trend: critical-slowing guard
TAU_CHUNK  = 200.0   # ~28 pitch periods at omega* ~ 1 (nondimensional tau)


def classify_orbit(rhs, args, y0, tau_chunk=TAU_CHUNK, floor_deg=FLOOR_DEG,
                   rel_tol=REL_TOL, max_chunks=MAX_CHUNKS):
    """March in chunks; decide LCO vs collapse.

    Returns (amp_deg | None, y_final, n_chunks, flag). amp is the late-window
    peak |pitch| in degrees; None means collapsed (or monotonically dying).
    """
    y = np.asarray(y0, float).copy()
    amps = []
    for k in range(max_chunks):
        sol = solve_ivp(rhs, (0.0, tau_chunk), y, args=args, method="LSODA",
                        rtol=1e-7, atol=1e-9,
                        t_eval=np.linspace(0.4 * tau_chunk, tau_chunk, 2000))
        y = sol.y[:, -1].copy()
        amps.append(float(np.degrees(np.max(np.abs(sol.y[1])))))
        if amps[-1] < floor_deg:
            return None, y, k + 1, "floor"
        if k >= 1 and abs(amps[-1] - amps[-2]) <= rel_tol * amps[-1]:
            return amps[-1], y, k + 1, "conv"
        if k >= 2 and amps[-1] < 3.0 and amps[-1] < amps[-2] < amps[-3]:
            return None, y, k + 1, "decay<3"    # small and monotone dying
    dec = (amps[-2] - amps[-1]) / max(amps[-1], 1e-9)
    if dec > 0.05:
        return None, y, max_chunks, "decay~"
    return amps[-1], y, max_chunks, "slow~"


def down_sweep(rhs, args_of_U, y_seed, U_hi, dU, U_min, refine,
               tau_chunk=TAU_CHUNK, verbose=True):
    """Warm-started down-sweep. Returns (branch, U_fold, amp_at_fold, seconds).

    branch is the list of sustained (U, amp) points. U_fold is the last speed
    with a sustained LCO, bisected to `refine` m/s. If the LCO survives to
    U_min the sweep stops there (deep branch; U_fold <= U_min).
    """
    t0 = time.perf_counter()
    amp, y, k, fl = classify_orbit(rhs, args_of_U(U_hi), y_seed, tau_chunk)
    if amp is None:
        raise RuntimeError("no LCO at sweep start -- raise U_hi or the seed")
    if verbose:
        print(f"    U={U_hi:7.3f}  amp={amp:6.2f}  [{fl},{k}ch]", flush=True)
    branch, U_last, y_last = [(U_hi, amp)], U_hi, y
    U = U_hi
    while U - dU >= U_min - 1e-9:
        U = round(U - dU, 6)
        amp, yf, k, fl = classify_orbit(rhs, args_of_U(U), y_last, tau_chunk)
        if verbose:
            shown = "DEAD" if amp is None else f"{amp:6.2f}"
            print(f"    U={U:7.3f}  amp={shown}  "
                  f"[{fl},{k}ch,{time.perf_counter()-t0:6.1f}s]", flush=True)
        if amp is None:
            lo, hi, y_hi, amp_hi = U, U_last, y_last, branch[-1][1]
            while hi - lo > refine:
                mid = round(0.5 * (lo + hi), 4)
                a2, y2, k2, f2 = classify_orbit(rhs, args_of_U(mid), y_hi,
                                                tau_chunk)
                if verbose:
                    shown = "DEAD" if a2 is None else f"{a2:6.2f}"
                    print(f"      bisect U={mid:7.3f}  amp={shown}  "
                          f"[{f2},{k2}ch]", flush=True)
                if a2 is None:
                    lo = mid
                else:
                    hi, y_hi, amp_hi = mid, y2, a2
                    branch.append((mid, a2))
            return branch, hi, amp_hi, time.perf_counter() - t0
        branch.append((U, amp))
        U_last, y_last = U, yf
    return branch, U_last, branch[-1][1], time.perf_counter() - t0