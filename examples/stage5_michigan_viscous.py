"""Stage 5b -- MICHIGAN wing, viscous rungs (QS+Stall, ONERA).

The static rung produces NOTHING: the wing trims at its polar's Cl_max
(local slope 0.63/rad vs 2pi -- a 90% collapse), so the flutter engine has
no fuel; the null is spot-checked here and fully tested in the stage-3
record. The ONERA rung, with Petot's published NACA 0012 laws and nothing
tuned, is the only rung on either axis that is subcritical: fold 9.594,
Hopf 10.875 (read from the committed stage-4 record; the 490 s sweep is not
re-run). The inviscid branches (5a) are drawn in grey for reference; the
experimental baseline sits behind everything.

Idempotent per rung (FORCE=1 recomputes the stall spot-check).
Writes: keys QS_Stall, ONERA in results/michigan_branches.json
        + results/F_michigan_viscous.png
Run:    PYTHONPATH=src python examples/stage5_michigan_viscous.py   (~15 s)
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
from mflco.aero.qs_stall import QSStall, load_polar
from stage5_common import make_classify, merge_write, COLORS
from exp_digitised import draw_experiment

OUT = "results/michigan_branches.json"
BETA = 1.326
ONERA_HOPF = 10.875

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
    d = rec["ONERA"]
    c = COLORS["ONERA"]
    pts = [t for t in sorted(map(tuple, d["branch"])) if t[1] <= 26.0]
    u, a = zip(*pts)
    ax.plot(u, a, "o-", ms=3.5, lw=1.5, color=c,
            label="ONERA, published laws (subcritical)")
    ax.plot([d["U_fold"]], [d["amp_at_fold"]], "v", ms=9, color=c)
    ax.annotate("display truncated at polar validity (~25 deg)",
                (u[-1], a[-1]), textcoords="offset points",
                xytext=(-152, 8), fontsize=6.8, color=c)
    ax.text(10.6, 28.0, "QS+Stall: equilibrium stable everywhere\n"
            "(no Hopf, no LCO -- tested null)",
            fontsize=7.3, color=COLORS["QS_Stall"])
    ax.set_title("Michigan wing, viscous axis: the static rung is silent; "
                 "only dynamic stall is subcritical", fontsize=9.5)
    ax.set_xlabel("airspeed U [m/s]")
    ax.set_ylabel("pitch LCO amplitude [deg]")
    ax.set_xlim(8.5, 14.6)
    ax.set_ylim(-1.5, 31)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7.5, loc="center left")
    plt.tight_layout()
    plt.savefig("results/F_michigan_viscous.png", dpi=160)


if __name__ == "__main__":
    rec = json.load(open(OUT)) if os.path.exists(OUT) else {}
    force = bool(os.environ.get("FORCE"))
    print("Michigan wing -- viscous axis", flush=True)
    upd = {}
    if force or "QS_Stall" not in rec:
        alpha_deg, cl = load_polar()
        aero = QSStall(p, 0.0, alpha_deg=alpha_deg, cl=cl, trim_deg=10.0)
        verdicts = {}
        for U in (10.5, 12.5, 13.5):
            y = np.zeros(4)
            y[1] = np.radians(10.0)
            a, _ = cls(aero, conv(U), y)
            verdicts[U] = None if a is None else round(a, 2)
        assert all(v is None for v in verdicts.values()), \
            f"stall rung sustained an orbit?! {verdicts}"
        upd["QS_Stall"] = dict(hopf=None, branch=[],
                               spot_checks=list(verdicts),
                               note="no oscillation at any probe; full "
                                    "tested null in the stage-3 record")
        print(f"  QS+Stall null at U = {list(verdicts)}", flush=True)
    if force or "ONERA" not in rec:
        src = "results/stage4_onera_published_laws.json"
        assert os.path.exists(src), "run examples/stage4_onera.py first"
        n = json.load(open(src))["nominal"]
        upd["ONERA"] = dict(hopf=ONERA_HOPF, branch=n["branch"],
                            U_fold=n["U_fold"],
                            amp_at_fold=n["amp_at_fold"],
                            note="from the committed stage-4 record")
        print(f"  ONERA fold {n['U_fold']}, Hopf {ONERA_HOPF} (from {src})",
              flush=True)
    rec = merge_write(OUT, upd)
    figure(rec)
    print(f"-> {OUT} + results/F_michigan_viscous.png", flush=True)