"""Shared machinery for the stage-5 branch runners (5a-5d).

One classify wrapper with runaway guards, the warm-started down-sweep (fold
by bracket-and-bisect), the up-sweep (supercritical branches), the small-kick
Hopf bisection with an auto-expanding upper bracket, the schematic 0015-like
polar used by the Bristol viscous rungs, and the merge-write for the per-rig
result JSONs (each runner owns its keys; reruns of one file never clobber
the other's record).
"""
import json
import os
import time

import numpy as np

from mflco.model.eom import structural_rhs
from mflco.model.fold_detector import classify_orbit


def make_classify(p):
    def _classify(aero, us, y, cap=45.0):
        if not np.all(np.isfinite(y)):
            return None, y
        try:
            amp, yf, _, _ = classify_orbit(structural_rhs, (p, aero, us), y)
        except Exception:
            return None, y
        if amp is not None and (amp > cap or not np.all(np.isfinite(yf))):
            return None, y
        return amp, (yf if amp is not None else y)
    return _classify


def down_sweep(cls, conv, name, aero, seeds, dU=0.25, U_min=13.5, refine=0.05):
    t0 = time.perf_counter()
    amp = None
    for U0, kick in seeds:
        y = np.zeros(4 + getattr(aero, "n_aero_states", 0))
        y[1] = np.radians(kick)
        amp, y = cls(aero, conv(U0), y)
        if amp is not None:
            break
    assert amp is not None, f"{name}: no seed caught"
    branch = [(U0, amp)]
    U = U0
    while U - dU >= U_min:
        U = round(U - dU, 3)
        a, y2 = cls(aero, conv(U), y)
        if a is None:
            lo, hi, yh = U, branch[-1][0], y
            while hi - lo > refine:
                mid = round(0.5 * (lo + hi), 3)
                a2, y2 = cls(aero, conv(mid), yh)
                if a2 is None:
                    lo = mid
                else:
                    hi, yh = mid, y2
                    branch.append((mid, a2))
            branch.sort()
            print(f"  {name:22s} fold {hi:6.2f}   amp_at_fold "
                  f"{branch[0][1]:5.2f} deg   "
                  f"[{time.perf_counter()-t0:4.0f} s]", flush=True)
            return dict(U_fold=hi,
                        branch=[(round(u, 2), round(a, 2))
                                for u, a in branch])
        branch.append((U, a))
        y = y2
    branch.sort()
    return dict(U_fold=None,
                branch=[(round(u, 2), round(a, 2)) for u, a in branch])


def up_sweep(cls, conv, name, aero, U0, kick, dU, U_max, amp_stop=30.0):
    t0 = time.perf_counter()
    y = np.zeros(4 + getattr(aero, "n_aero_states", 0))
    y[1] = np.radians(kick)
    amp, y = cls(aero, conv(U0), y)
    assert amp is not None, f"{name}: seed dead at {U0}"
    br = [(U0, amp)]
    U = U0
    while amp is not None and amp < amp_stop and U + dU <= U_max:
        U = round(U + dU, 3)
        amp, y2 = cls(aero, conv(U), y)
        if amp is None:
            break
        br.append((U, amp))
        y = y2
    print(f"  {name:22s} {br[0][0]:.2f}({br[0][1]:.1f} deg) -> "
          f"{br[-1][0]:.2f}({br[-1][1]:.1f} deg)   "
          f"[{time.perf_counter()-t0:4.0f} s]", flush=True)
    return dict(U_fold=None,
                branch=[(round(u, 2), round(a, 2)) for u, a in br])


def hopf_bisect(cls, conv, name, aero, lo, hi, kick=0.4, tol=0.1):
    """Small-kick stability: decay = stable equilibrium; capture = unstable."""
    t0 = time.perf_counter()

    def unstable(U):
        y = np.zeros(4 + getattr(aero, "n_aero_states", 0))
        y[1] = np.radians(kick)
        a, _ = cls(aero, conv(U), y)
        return a is not None

    assert not unstable(lo), f"{name}: lower bracket already unstable"
    while not unstable(hi):
        hi = round(hi + 1.0, 3)
        assert hi <= 30.0, f"{name}: no instability below 30 m/s"
    while hi - lo > tol:
        mid = round(0.5 * (lo + hi), 3)
        lo, hi = (lo, mid) if unstable(mid) else (mid, hi)
    print(f"  {name:22s} Hopf {0.5*(lo+hi):6.2f}  (bracket {lo}-{hi})"
          f"   [{time.perf_counter()-t0:4.0f} s]", flush=True)
    return round(0.5 * (lo + hi), 2)


def bristol_schematic_polar():
    """Schematic 0015-like polar: slope 6.0/rad, stall onset 11 deg.
    Labelled schematic wherever used; a digitised table is future work."""
    ag = np.array([-180, -160, -120, -90, -60, -30, -20, -14, -12, -11, -8,
                   -4, 0, 4, 8, 11, 12, 14, 20, 30, 60, 90, 120, 160, 180],
                  float)

    def cl15(a):
        s = np.sign(a)
        x = abs(a)
        if x <= 11:
            return 6.0 * np.radians(a)
        if x <= 14:
            return s * (1.10 - 0.02 * (x - 11))
        if x <= 20:
            return s * (1.04 - 0.03 * (x - 14))
        if x <= 90:
            return s * (0.45 + 0.225 * np.sin(np.radians(2 * x))
                        / np.sin(np.radians(40)))
        return s * 0.1 * np.sin(np.radians(2 * (180 - x)))

    return ag, np.array([cl15(a) for a in ag])


def merge_write(path, updates):
    """Update this runner's keys in the per-rig JSON without clobbering."""
    rec = json.load(open(path)) if os.path.exists(path) else {}
    rec.update(updates)
    json.dump(rec, open(path, "w"), indent=1)
    return rec


def zero_line(ax, x_lo, x_hi, fold, hopf, c):
    """Three-regime equilibrium line at y = 0 (overlap is fine)."""
    if hopf is None:
        ax.plot([x_lo, x_hi], [0, 0], "-", lw=1.1, color=c, alpha=0.75)
        return
    knee = fold if fold is not None else hopf
    ax.plot([x_lo, knee], [0, 0], "-", lw=1.1, color=c, alpha=0.75)
    if fold is not None and fold < hopf:
        ax.plot([fold, hopf], [0, 0], ":", lw=1.5, color=c, alpha=0.9)
    ax.plot([hopf, x_hi], [0, 0], "--", lw=1.0, color=c, alpha=0.65)
    ax.plot([hopf], [0], "o", ms=6, color=c, mfc="none", mec=c, mew=1.5)


COLORS = {"QS": "tab:blue", "QS_Stall": "tab:green",
          "Peters": "tab:orange", "ONERA": "tab:purple"}