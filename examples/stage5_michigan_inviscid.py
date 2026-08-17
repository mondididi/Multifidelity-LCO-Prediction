"""Stage 5a -- MICHIGAN wing, inviscid rungs (QS, Peters N=6).

Both supercritical: no fold exists on this axis. QS flutters at 9.665
(-26.7%: no wake lag, pitch over-destabilised) with a near-vertical branch;
Peters at 13.148 (-0.32%, UNFITTED) with the beta-anchored branch through
14 deg at 13.5. The experimental baseline (digitised, pitch) is drawn
behind: its fold at 11.85 is exactly where neither rung can go.

Beware: section_from_params() gives beta = 0; the cubic is set explicitly
(beta = 1.326, the stage-2 anchor), as in examples/time_costs.py.

Idempotent per rung (FORCE=1 recomputes).
Writes: keys QS, Peters in results/michigan_branches.json
        + results/F_michigan_inviscid.png
Run:    PYTHONPATH=src python examples/stage5_michigan_inviscid.py  (~1 min)
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from mflco.model.params import TypicalSectionParameters
from mflco.model.michigan_params import calibrate_michigan, structural_zeta
from mflco.aero.quasi_steady import QuasiSteady
from mflco.aero.peters_finite import PetersFinite
from stage5_common import (make_classify, up_sweep, merge_write,
                           COLORS)
from exp_digitised import draw_experiment

OUT = "results/michigan_branches.json"
BETA = 1.326
HOPF = {"QS": 9.665, "Peters": 13.148}     # linear flutter speeds (locked)

cal = calibrate_michigan(zeta=structural_zeta())
p = TypicalSectionParameters(
    a=cal.a, x_alpha=cal.x_alpha, r_alpha_sq=cal.r_alpha_sq,
    omega_ratio=cal.omega_ratio, mu=cal.mu, beta=BETA, delta=0.0,
    zeta_h=cal.zeta, zeta_alpha=cal.zeta)
conv = lambda U: float(cal.ms_to_ustar(U))
cls = make_classify(p)


def figure(rec):
    fig, ax = plt.subplots(figsize=(7.2, 4.3))
    draw_experiment(ax, "michigan", x_max=14.6, label_prefix="experiment")
    for key, lab in (("QS", "QS (supercritical)"),
                     ("Peters", "Peters N=6 (supercritical)")):
        c = COLORS[key]
        u, a = zip(*sorted(map(tuple, rec[key]["branch"])))
        ax.plot(u, a, "o-", ms=3.5, lw=1.5, color=c, label=lab)
    ax.set_title("Michigan wing, inviscid axis: both rungs supercritical -- "
                 "no fold exists here", fontsize=9.5)
    ax.set_xlabel("airspeed U [m/s]")
    ax.set_ylabel("pitch LCO amplitude [deg]")
    ax.set_xlim(8.5, 14.6)
    ax.set_ylim(-1.5, 31)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7.5, loc="upper left")
    plt.tight_layout()
    plt.savefig("results/F_michigan_inviscid.png", dpi=160)


if __name__ == "__main__":
    rec = json.load(open(OUT)) if os.path.exists(OUT) else {}
    force = bool(os.environ.get("FORCE"))
    print("Michigan wing -- inviscid axis", flush=True)
    upd = {"beta": BETA}
    if force or "QS" not in rec:
        d = up_sweep(cls, conv, "QS(super)", QuasiSteady(p, 0.0),
                     9.8, 1.5, 0.1, 12.0)
        d["hopf"] = HOPF["QS"]
        d["note"] = "supercritical; truncated ~30 deg (model validity)"
        upd["QS"] = d
    if force or "Peters" not in rec:
        aero = PetersFinite(p, 0.0, N=6)
        y0 = np.zeros(4 + aero.n_aero_states)
        y0[1] = np.radians(8.0)
        a0, y0 = cls(aero, conv(13.5), y0)
        assert a0 and 12.0 < a0 < 16.0, f"anchor check failed: {a0}"
        print(f"  Peters anchor 13.5 -> {a0:.2f} deg (expect ~14)",
              flush=True)
        br = [(13.5, a0)]
        U, y = 13.5, y0.copy()
        while U - 0.1 >= 13.0:
            U = round(U - 0.1, 3)
            a, y2 = cls(aero, conv(U), y)
            if a is None:
                break
            br.append((U, a))
            y = y2
        U, y = 13.5, y0.copy()
        while U + 0.25 <= 14.25:
            U = round(U + 0.25, 3)
            a, y2 = cls(aero, conv(U), y)
            if a is None:
                break
            br.append((U, a))
            y = y2
        br.sort()
        upd["Peters"] = dict(
            hopf=HOPF["Peters"],
            branch=[(round(u, 2), round(a, 2)) for u, a in br],
            note="supercritical from 13.148 (flutter, unfitted)")
        print(f"  Peters {br[0][0]}({br[0][1]:.1f}) -> "
              f"{br[-1][0]}({br[-1][1]:.1f})", flush=True)
    rec = merge_write(OUT, upd)
    figure(rec)
    print(f"-> {OUT} + results/F_michigan_inviscid.png", flush=True)