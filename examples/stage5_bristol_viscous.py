"""Stage 5d -- BRISTOL wing, viscous rungs (QS+Stall, ONERA).

The test of the structural verdict: if the Bristol fold belongs to the
spring, separation physics should change little -- and it is nearly inert
by construction, because at zero trim on a symmetric section the stall
deficit is <0.05 below the ~11 deg onset (vs 0.62 at 15 deg). Both viscous
rungs use the SCHEMATIC 0015-like polar (slope 6.0/rad, stall 11 deg;
digitised table = future work). Record: QS+Stall fold 14.94 (+0.25 vs QS),
ONERA fold 17.34 (+0.43 vs Peters); Hopfs 24.07 and 24.92. Inviscid
branches (5c) in grey for the shift; experiment behind on the twin axis.

Idempotent per rung (FORCE=1 recomputes).
Writes: keys QS_Stall, ONERA in results/bristol_branches.json
        + results/F_bristol_viscous.png
Run:    PYTHONPATH=src python examples/stage5_bristol_viscous.py  (~4 min)
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from mflco.model.barton_params import BartonCal
from mflco.aero.qs_stall import QSStall
from mflco.aero.onera_stall import ONERAStall
from stage5_common import (make_classify, down_sweep, hopf_bisect,
                           merge_write, bristol_schematic_polar,
                           COLORS)
from exp_digitised import draw_experiment

OUT = "results/bristol_branches.json"

cal = BartonCal(2)
p = cal.section()
conv = lambda U: float(cal.ms_to_ustar(U))
cls = make_classify(p)
ALPHA, CL = bristol_schematic_polar()


def figure(rec):
    fig, ax = plt.subplots(figsize=(7.2, 4.3))
    ax2 = ax.twinx()
    draw_experiment(ax2, "bristol", x_max=26.0,
                    label_prefix="experiment, heave [mm]")
    ax2.set_ylim(-6, 132)
    ax2.set_ylabel("experimental heave [mm] (right axis)", fontsize=8.5)
    for key, lab in (("QS_Stall", "QS+Stall (schematic polar)"),
                     ("ONERA", "ONERA (schematic polar)")):
        c = COLORS[key]
        u, a = zip(*sorted(map(tuple, rec[key]["branch"])))
        ax.plot(u, a, "o-", ms=3.5, lw=1.5, color=c, label=lab)
        ax.plot([u[0]], [a[0]], "v", ms=8, color=c)
    ax.set_title("Bristol wing, viscous axis: separation physics is nearly "
                 "inert (+0.25 / +0.43 m/s)", fontsize=9.5)
    ax.set_xlabel("airspeed U [m/s]")
    ax.set_ylabel("model pitch LCO amplitude [deg]")
    ax.set_xlim(13, 26)
    ax.set_ylim(-1.4, 31)
    ax.grid(alpha=0.25)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=7.2, loc="upper left")
    plt.tight_layout()
    plt.savefig("results/F_bristol_viscous.png", dpi=160)


if __name__ == "__main__":
    rec = json.load(open(OUT)) if os.path.exists(OUT) else {}
    force = bool(os.environ.get("FORCE"))
    print("Bristol wing -- viscous axis", flush=True)
    upd = {"polar": "schematic 0015-like (slope 6.0/rad, stall 11 deg)"}
    if force or "QS_Stall" not in rec:
        st = QSStall(p, 0.0, alpha_deg=ALPHA, cl=CL, trim_deg=0.0)
        d = down_sweep(cls, conv, "QS+Stall(schematic)", st,
                       [(17.0, 9.0), (18.5, 11.0)])
        d["hopf"] = hopf_bisect(cls, conv, "QS+Stall hopf", st, 15.5, 20.0)
        upd["QS_Stall"] = d
    if force or "ONERA" not in rec:
        on = ONERAStall(p, 0.0, N=6, alpha_deg=ALPHA, cl=CL, trim_deg=0.0)
        d = down_sweep(cls, conv, "ONERA(schematic)", on, [(22.0, 15.0)],
                       U_min=14.0)
        d["hopf"] = hopf_bisect(cls, conv, "ONERA hopf", on, 22.0, 26.5)
        upd["ONERA"] = d
    rec = merge_write(OUT, upd)
    figure(rec)
    print(f"-> {OUT} + results/F_bristol_viscous.png", flush=True)