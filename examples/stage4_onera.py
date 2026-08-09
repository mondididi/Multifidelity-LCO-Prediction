"""Stage 4 -- ONERA dynamic stall: does dynamic separation produce the fold?

The hinge run of the viscous axis. ONERAStall = Peters (attached, validated
-0.32% on flutter) + the ONERA second-order stalled-load equation, so the ONLY
physics added relative to the Peters rung is dynamic separation with
hysteresis. Structure held FIXED at the stage-2 build (delta = 0, beta = 1.326)
-- the ablation discipline of stage 3.

Stability is characterised EMPIRICALLY (kick probes + the fold detector)
rather than by eigenvalue sweep: the stalled-load states are not in the
descriptor matrices, and the deficit slope at trim (dDelta/dalpha ~ 5.65 /rad,
90% of 2*pi) makes the stall equation linearly active, so the honest
linearisation is the time-marching one. Per airspeed:

    small kick (0.5 deg) grows  -> linearly unstable there (Hopf below)
    small dead + large kick LCO -> bistable there (subcritical evidence)
    both dead                   -> stable

If any LCO is found, the fold detector reads the fold by down-sweep, exactly
as for every other rung.

COEFFICIENT ROBUSTNESS ("sweep" phase): a_stall/r_stall/e_stall are borrowed
NACA 0012 constants standing in for the paper's Delta-dependent splines, so no
single triple is trusted. A 3x3x3 grid spanning the plausible ranges is
screened at three airspeeds; the conclusion of this stage is the SET of
outcomes, not one run. Pre-registered rule (before any result was seen): a
fold produced is strong evidence dynamic separation suffices; a fold not
produced is weak evidence, because borrowed coefficients cannot distinguish
"dynamic stall is insufficient" from "wrong numbers for this airfoil".

Cold run: nominal ~2 min, sweep ~5 min; cached incrementally to the system
temp dir -- rerun is instant, delete the JSON to recompute.

Run:  PYTHONPATH=src python examples/stage4_onera.py nominal|sweep|report
"""
import itertools
import json
import os
import sys
import tempfile
import time

import numpy as np

PHASE = sys.argv[1] if len(sys.argv) > 1 else "nominal"

from mflco.model.params import TypicalSectionParameters
from mflco.model.michigan_params import calibrate_michigan, structural_zeta
from mflco.model.eom import structural_rhs
from mflco.model.fold_detector import classify_orbit, down_sweep
from mflco.aero.onera_stall import ONERAStall, A_STALL, R_STALL, E_STALL

# --- constants -------------------------------------------------------------------
RIG_FLUTTER = 13.19     # Garcia Perez J063736
RIG_FOLD    = 11.85     # Garcia Perez turning point
N_INFLOW    = 6         # Peters states, as stages 1-2
BETA_S2     = 1.326     # stage-2 structural build, held fixed (delta = 0)
KICK_SMALL  = 0.5       # deg -- linear-stability probe (floor is also 0.5)
KICK_LARGE  = [8.0, 15.0]      # deg -- basin probes for finite-amplitude LCOs
U_NOMINAL   = [6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 11.85, 12.5, 13.19,
               14.0, 15.0, 16.0]
U_SWEEP     = [8.0, 11.85, 16.0]               # low island / rig fold region / upper
GRID = dict(a_stall=[0.10, 0.50],              # stalled-mode damping (corners)
            r_stall=[0.10, 0.40],              # stalled-mode stiffness (corners)
            e_stall=[0.25, 1.00])              # deficit-rate gain (corners)
CACHE = os.path.join(tempfile.gettempdir(), "stage4_onera.json")

cal = calibrate_michigan(zeta=structural_zeta())


def _section():
    return TypicalSectionParameters(
        a=cal.a, x_alpha=cal.x_alpha, r_alpha_sq=cal.r_alpha_sq,
        omega_ratio=cal.omega_ratio, mu=cal.mu, beta=BETA_S2, delta=0.0,
        zeta_h=cal.zeta, zeta_alpha=cal.zeta)


def _aero(p, **kw):
    return ONERAStall(p, 0.0, N=N_INFLOW, **kw)


def _probe(p, aero, U_ms, kick_deg):
    """classify_orbit from a pure-pitch kick.

    Returns (amp_deg | None | "BLOWUP", flag). BLOWUP marks an orbit that left
    the finite domain -- unbounded response at that corner, recorded as its own
    outcome rather than crashing the sweep (seen at low-damping/high-gain
    corners where the stalled mode pumps beyond the polar's validity)."""
    y0 = np.zeros(4 + aero.n_aero_states)
    y0[1] = np.radians(kick_deg)
    try:
        amp, _, k, fl = classify_orbit(
            structural_rhs, (p, aero, float(cal.ms_to_ustar(U_ms))), y0)
    except (ValueError, FloatingPointError):
        return "BLOWUP", "blowup"
    return amp, f"{fl},{k}ch"


def _characterize(p, aero, U_list, verbose=True):
    """Per-U outcomes: small-kick growth and large-kick capture."""
    rows = []
    for U in U_list:
        s_amp, s_fl = _probe(p, aero, U, KICK_SMALL)
        l_amp, l_fl = None, "-"
        if s_amp is None:                       # (BLOWUP skips basin hunt too)                       # only hunt basins if eq. stable
            kicks = KICK_LARGE if PHASE == "nominal" else KICK_LARGE[:1]
            for kk in kicks:                    # sweep phase: 8 deg only (cost)
                l_amp, l_fl = _probe(p, aero, U, kk)
                if l_amp is not None:
                    break
        rows.append(dict(U=U, small=s_amp, large=l_amp))
        if verbose:
            s = ("BLOWUP" if s_amp == "BLOWUP"
                 else "grows" if s_amp is not None else "dead")
            l = f"LCO {l_amp:.1f}" if l_amp is not None else \
                ("-" if s_amp is not None else "dead")
            print(f"    U={U:6.2f}  small({KICK_SMALL}): {s:5s}  "
                  f"large: {l}", flush=True)
    return rows


def _load():
    try:
        return json.load(open(CACHE))
    except (FileNotFoundError, json.JSONDecodeError):
        return {}


def _save(c):
    json.dump(c, open(CACHE, "w"), indent=1)


# --- phases ----------------------------------------------------------------------
cache = _load()

if PHASE == "nominal":
    print(f"== ONERA nominal (a={A_STALL}, r={R_STALL}, e={E_STALL}) ==",
          flush=True)
    p = _section()
    aero = _aero(p)
    nom = cache.setdefault("nominal", {"rows": [], "probe_seconds": 0.0})
    done = {r["U"] for r in nom["rows"]}
    for U in U_NOMINAL:
        if U in done:
            continue
        t0 = time.perf_counter()
        row = _characterize(p, aero, [U])[0]
        nom["rows"].append(row)
        nom["probe_seconds"] = round(nom["probe_seconds"]
                                     + time.perf_counter() - t0, 1)
        _save(cache)
    rows = nom["rows"]
    lco_at = [r for r in rows if (r["small"] or r["large"])]
    if lco_at and "U_fold" not in nom:
        # resumable warm-started down-sweep: same protocol as
        # fold_detector.down_sweep, checkpointing (U, amp, state) per step so
        # an interrupted run continues from the last converged orbit
        argU = lambda U: (p, aero, float(cal.ms_to_ustar(U)))
        sw = nom.setdefault("sweep_state", {})
        if not sw:
            U_top = max(r["U"] for r in lco_at)
            y0 = np.zeros(4 + aero.n_aero_states)
            y0[1] = np.radians(8.0)
            print(f"\ndown-sweep (resumable) from U={U_top}", flush=True)
            amp, y, k, fl = classify_orbit(structural_rhs, argU(U_top), y0)
            assert amp is not None, "seed failed to capture the LCO"
            sw.update(U=U_top, amp=round(amp, 2), y=list(map(float, y)),
                      branch=[(U_top, round(amp, 2))], secs=0.0)
            _save(cache)
        t0 = time.perf_counter()
        while True:
            U_next = round(sw["U"] - 0.25, 3)
            if U_next < 4.0:
                U_dead = None
                break
            amp, y, k, fl = classify_orbit(structural_rhs, argU(U_next),
                                           np.asarray(sw["y"]))
            sw["secs"] = round(sw["secs"] + time.perf_counter() - t0, 1)
            t0 = time.perf_counter()
            print(f"    U={U_next:7.3f}  "
                  f"amp={'DEAD' if amp is None else f'{amp:6.2f}'}  "
                  f"[{fl},{k}ch, {sw['secs']:.0f}s]", flush=True)
            if amp is None:
                U_dead = U_next
                break
            sw.update(U=U_next, amp=round(amp, 2), y=list(map(float, y)))
            sw["branch"].append((U_next, round(amp, 2)))
            _save(cache)
        if U_dead is not None:
            lo, hi = U_dead, sw["U"]                 # bisect to 0.05 m/s
            y_hi = np.asarray(sw["y"])
            while hi - lo > 0.05:
                mid = 0.5 * (lo + hi)
                amp, y, k, fl = classify_orbit(structural_rhs, argU(mid), y_hi)
                print(f"      bisect U={mid:7.3f}  "
                      f"amp={'DEAD' if amp is None else f'{amp:6.2f}'}",
                      flush=True)
                if amp is None:
                    lo = mid
                else:
                    hi, y_hi = mid, y
                    sw.update(U=mid, amp=round(amp, 2))
            nom.update(U_fold=round(hi, 3), amp_at_fold=sw["amp"],
                       sweep_seconds=sw["secs"],
                       branch=sorted(sw["branch"]))
            _save(cache)
            print(json.dumps({k: v for k, v in nom.items()
                              if k not in ("rows", "branch", "sweep_state")},
                             indent=1))
    elif not lco_at:
        print("\nno self-sustained LCO at any probed speed/kick "
              "(nominal coefficients)")

elif PHASE == "island":
    # the low-speed instability island (unstable 7-9, bistable at 6): its fold
    # is the LOWEST speed with sustained oscillation, i.e. the constraint this
    # model actually reports. Seed on the island, warm-sweep down.
    print("== ONERA island fold (down-sweep from 8.0) ==", flush=True)
    p = _section()
    aero = _aero(p)
    y0 = np.zeros(4 + aero.n_aero_states)
    y0[1] = np.radians(KICK_SMALL)              # grows onto the island LCO
    argU = lambda U: (p, aero, float(cal.ms_to_ustar(U)))
    br, U_fold, amp_f, secs = down_sweep(
        structural_rhs, argU, y0, 8.0, 0.25, 3.0, 0.05)
    cache["island"] = dict(U_fold=round(U_fold, 3),
                           amp_at_fold=round(amp_f, 2),
                           sweep_seconds=round(secs, 1),
                           branch=[(round(u, 3), round(a, 2)) for u, a in br])
    _save(cache)
    print(json.dumps({k: v for k, v in cache["island"].items()
                      if k != "branch"}, indent=1))

elif PHASE == "sweep":
    print("== ONERA coefficient robustness sweep (3x3x3) ==", flush=True)
    p = _section()
    sw = cache.setdefault("sweep", {})
    for a_s, r_s, e_s in itertools.product(*GRID.values()):
        key = f"a{a_s}_r{r_s}_e{e_s}"
        if key in sw:
            continue
        aero = _aero(p, a_stall=a_s, r_stall=r_s, e_stall=e_s)
        rows = _characterize(p, aero, U_SWEEP, verbose=False)
        outcomes = [v for r in rows for v in (r["small"], r["large"])]
        verdict = ("BLOWUP" if "BLOWUP" in outcomes
                   else "LCO" if any(isinstance(v, float) for v in outcomes)
                   else "none")
        sw[key] = dict(rows=rows, verdict=verdict)
        _save(cache)
        print(f"  {key:22s} -> {verdict}", flush=True)

elif PHASE == "report":
    nom = cache.get("nominal", {})
    print("NOMINAL:", "fold at " + str(nom.get("U_fold")) if "U_fold" in nom
          else "no LCO found")
    sw = cache.get("sweep", {})
    n_lco = sum(1 for v in sw.values() if v["verdict"] == "LCO")
    print(f"SWEEP: {n_lco}/{len(sw)} coefficient triples produce any LCO")
    for k, v in sw.items():
        if v["verdict"] == "LCO":
            print(f"  {k}: "
                  + "; ".join(f"U={r['U']}: small="
                              f"{'grow' if r['small'] else 'dead'}, large="
                              f"{r['large'] or '-'}" for r in v["rows"]))
    print(f"cache -> {CACHE}")
