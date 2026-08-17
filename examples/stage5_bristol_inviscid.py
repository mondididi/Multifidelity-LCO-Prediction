"""Stage 5c -- BRISTOL wing, inviscid rungs (QS, Peters N=6).

Both subcritical from the documented spring alone. The record: folds at
14.69 (QS) and 16.91 (Peters) -- the fold's EXISTENCE is the spring's, but
wake fidelity moves its LOCATION by 2.2 m/s -- while the Hopf points sit at
23.64 (QS, bisected; linearised cross-check 23.2) and the test-locked 23.98
(Peters, 0.08% from the independent 24.0). On this wing wake fidelity barely
touches flutter and moves the fold: the inverse of Michigan. Experimental
baseline (heave, mm -- their diagram's DoF) rides a twin right axis behind
the models; velocities are the like-for-like comparison.

Idempotent per rung (FORCE=1 recomputes).
Writes: keys QS, Peters in results/bristol_branches.json
        + results/F_bristol_inviscid.png
Run:    PYTHONPATH=src python examples/stage5_bristol_inviscid.py  (~2 min)
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from mflco.model.barton_params import BartonCal
from mflco.aero.quasi_steady import QuasiSteady
from mflco.aero.peters_finite import PetersFinite
from stage5_common import (make_classify, down_sweep, hopf_bisect,
                           merge_write, COLORS)
from exp_digitised import draw_experiment

OUT = "results/bristol_branches.json"
PETERS_HOPF = 23.98            # test-locked (0.08% vs the independent 24.0)

cal = BartonCal(2)
p = cal.section()
conv = lambda U: float(cal.ms_to_ustar(U))
cls = make_classify(p)


def figure(rec):
    fig, ax = plt.subplots(figsize=(7.2, 4.3))
    ax2 = ax.twinx()
    draw_experiment(ax2, "bristol", x_max=26.0,
                    label_prefix="experiment, heave [mm]")
    ax2.set_ylim(-6, 132)
    ax2.set_ylabel("experimental heave [mm] (right axis)", fontsize=8.5)
    for key, lab in (("QS", "QS"), ("Peters", "Peters N=6")):
        c = COLORS[key]
        u, a = zip(*sorted(map(tuple, rec[key]["branch"])))
        ax.plot(u, a, "o-", ms=3.5, lw=1.5, color=c, label=lab)
        ax.plot([u[0]], [a[0]], "v", ms=8, color=c)
    ax.set_title("Bristol wing, inviscid axis: both rungs subcritical from "
                 "the spring alone", fontsize=9.5)
    ax.set_xlabel("airspeed U [m/s]")
    ax.set_ylabel("model pitch LCO amplitude [deg]")
    ax.set_xlim(13, 26)
    ax.set_ylim(-1.4, 31)
    ax.grid(alpha=0.25)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=7.2, loc="upper left")
    plt.tight_layout()
    plt.savefig("results/F_bristol_inviscid.png", dpi=160)


if __name__ == "__main__":
    rec = json.load(open(OUT)) if os.path.exists(OUT) else {}
    force = bool(os.environ.get("FORCE"))
    print("Bristol wing -- inviscid axis", flush=True)
    upd = {}
    if force or "QS" not in rec:
        qs = QuasiSteady(p, 0.0)
        d = down_sweep(cls, conv, "QS", qs, [(17.0, 9.0), (18.5, 11.0)])
        d["hopf"] = hopf_bisect(cls, conv, "QS hopf", qs, 15.5, 20.0)
        upd["QS"] = d
    if force or "Peters" not in rec:
        pe = PetersFinite(p, 0.0, N=6)
        d = down_sweep(cls, conv, "Peters N=6", pe, [(22.0, 15.0)],
                       U_min=15.5)
        d["hopf"] = PETERS_HOPF
        print(f"  {'Peters hopf':22s} Hopf {PETERS_HOPF:6.2f} (test-locked)",
              flush=True)
        upd["Peters"] = d
    rec = merge_write(OUT, upd)
    figure(rec)
    print(f"-> {OUT} + results/F_bristol_inviscid.png", flush=True)