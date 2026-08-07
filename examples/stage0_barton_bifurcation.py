"""Stage 0 -- Barton case-2 bifurcation diagram (the verification figure).

Self-contained: computes everything it needs -- linear flutter (eigen sweep),
stable LCO branches (fold_detector.down_sweep), and unstable-branch estimates
(kick-threshold bisection) -- for both aero flavours, then renders the
diagram. Results are cached incrementally to a JSON in the system temp dir,
so an interrupted run resumes and a rerun plots instantly; delete the JSON to
recompute from scratch. Cold run ~5 min (Peters dominates).

Unstable branch method: at each speed between fold and Hopf, bisect the
initial pitch kick between "decays to equilibrium" and "captured by the
stable LCO". The threshold approximates the unstable orbit amplitude -- the
numerical analogue of the experimental protocols (Tartaruga's CBC stabilises
the unstable branch by control; Garcia Perez applies perturbations of
increasing amplitude). Caveat for the caption: a pure-pitch kick threshold is
a basin-boundary estimate, not the exact unstable orbit.

Experimental overlays are SPEEDS only (Hopf ~24, fold ~16 m/s): the rig's
published amplitude axis is heave [mm], this figure's is pitch [deg].

"""
import json
import os
import tempfile

import numpy as np

from mflco.model.barton_params import BartonCal
from mflco.model.analysis import modal_analysis
from mflco.model.eom import structural_rhs
from mflco.model.fold_detector import classify_orbit, down_sweep
from mflco.aero.quasi_steady import QuasiSteady
from mflco.aero.peters_finite import PetersFinite

N_INFLOW  = 6                                     # Peters states, as stage 0
UNSTABLE_SPEEDS = {"qs": [16.0, 18.0, 20.0, 22.0, 23.0],
                   "peters": [17.5, 19.0, 21.0, 23.0]}
CACHE = os.path.join(tempfile.gettempdir(), "stage0_barton_bifurcation.json")
PNG   = os.path.join(tempfile.gettempdir(), "stage0_barton_bifurcation.png")

cal = BartonCal(2)                                # the published campaign


# --- helpers ---------------------------------------------------------------------
def _make_aero(kind, p):
    if kind == "peters":
        return PetersFinite(p, 0.0, N=N_INFLOW), N_INFLOW
    return QuasiSteady(p, 0.0), 0


def _flutter_ms(aero, us_lo=6.0, us_hi=10.0, n=300):
    """Linear flutter [m/s]: min-damping zero crossing on a moderate bracket."""
    p = aero.params
    prev, Uf = None, None
    for u in np.linspace(us_lo, us_hi, n):
        _, d = modal_analysis(p, 0.0, u, aero)
        m = float(np.min(d))
        if prev is not None and prev > 0.0 >= m:
            Uf = u
            break
        prev = m
    assert Uf is not None, "no flutter crossing in bracket"
    return float(cal.ustar_to_ms(Uf))


def _threshold(kind, U_ms, hi_deg, lo_deg=0.3, tol_deg=0.3):
    """Bisect the pure-pitch kick between decay and capture by the stable LCO."""
    p = cal.section()
    aero, n_a = _make_aero(kind, p)
    args = (p, aero, float(cal.ms_to_ustar(U_ms)))
    lo, hi = lo_deg, hi_deg
    while hi - lo > tol_deg:
        mid = 0.5 * (lo + hi)
        y0 = np.zeros(4 + n_a)
        y0[1] = np.radians(mid)
        amp, _, k, fl = classify_orbit(structural_rhs, args, y0)
        print(f"      U={U_ms:5.2f} kick={mid:6.2f} -> "
              f"{'LCO' if amp is not None else 'decay'} [{fl},{k}ch]", flush=True)
        if amp is None:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _load_cache():
    try:
        return json.load(open(CACHE))
    except (FileNotFoundError, json.JSONDecodeError):
        return {}


def _save_cache(c):
    json.dump(c, open(CACHE, "w"), indent=1)


# --- compute (load-or-run, incremental) ------------------------------------------
cache = _load_cache()
for kind in ("qs", "peters"):
    c = cache.setdefault(kind, {})
    if "hopf" not in c:
        p = cal.section()
        aero, _ = _make_aero(kind, p)
        c["hopf"] = round(_flutter_ms(aero), 3)
        print(f"[{kind}] linear flutter: {c['hopf']} m/s", flush=True)
        _save_cache(cache)
    if "branch" not in c:
        print(f"[{kind}] stable branch (down-sweep from "
              f"{c['hopf'] + 0.5:.2f} m/s)...", flush=True)
        p = cal.section()
        aero, n_a = _make_aero(kind, p)
        y0 = np.zeros(4 + n_a)
        y0[1] = np.radians(11.0)
        argU = lambda U: (p, aero, float(cal.ms_to_ustar(U)))
        br, U_fold, amp_f, secs = down_sweep(
            structural_rhs, argU, y0, round(c["hopf"] + 0.5, 2), 0.25,
            max(c["hopf"] - 10.0, 5.0), 0.05)
        c["branch"] = [(round(u, 3), round(a, 2)) for u, a in br]
        c["fold"], c["amp_at_fold"] = round(U_fold, 3), round(amp_f, 2)
        print(f"[{kind}] fold {c['fold']} m/s (amp {c['amp_at_fold']} deg) "
              f"[{secs:.0f}s]", flush=True)
        _save_cache(cache)
    if "unstable" not in c:
        print(f"[{kind}] unstable-branch kick thresholds...", flush=True)
        stable = dict(c["branch"])
        Ugrid = sorted(stable)
        pts = []
        for U in UNSTABLE_SPEEDS[kind]:
            hi = float(np.interp(U, Ugrid, [stable[u] for u in Ugrid])) - 0.5
            pts.append((U, round(_threshold(kind, U, hi_deg=hi), 2)))
            c["unstable"] = pts
            _save_cache(cache)
        print(f"[{kind}] thresholds: {c['unstable']}", flush=True)

# --- plot ------------------------------------------------------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(8.6, 5.4))
colors = {"qs": "tab:blue", "peters": "tab:orange"}
for kind in ("qs", "peters"):
    c, col = cache[kind], colors[kind]
    br = np.array(sorted(c["branch"]))
    ax.plot(br[:, 0], br[:, 1], "o-", ms=3.5, color=col,
            label=f"mflco {kind}: stable LCO "
                  f"(Hopf {c['hopf']:.2f}, fold {c['fold']:.2f})")
    ax.plot([c["fold"]], [c["amp_at_fold"]], "v", color=col, ms=9)
    ax.plot([c["hopf"]], [0.0], "*", color=col, ms=12)
    un = np.array(c["unstable"])
    seq = np.vstack([[c["fold"], c["amp_at_fold"]], un, [c["hopf"], 0.0]])
    ax.plot(seq[:, 0], seq[:, 1], "s--", ms=5, mfc="none", color=col,
            alpha=0.8, label=f"mflco {kind}: unstable (kick threshold)")
ax.axvline(24.0, color="r", ls="--", lw=1.2, label="Tartaruga rig Hopf ~24")
ax.axvline(16.0, color="r", ls=":", lw=1.4, label="Tartaruga rig fold ~16")
ax.set_xlabel("airspeed U [m/s]")
ax.set_ylabel("pitch LCO amplitude [deg]")
ax.set_title("Stage 0: Barton case 2 through the mflco pipeline -- "
             "subcritical bifurcation diagram")
ax.legend(fontsize=8, loc="upper left")
ax.grid(alpha=0.25)
ax.set_xlim(13.5, 26.0)
ax.set_ylim(-0.4, 13.0)
plt.tight_layout()
plt.savefig(PNG, dpi=170)
print(f"\ncache  -> {CACHE}")
print(f"figure -> {PNG}")