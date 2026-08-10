"""Build the paper's figure set and results matrix into results/.

Regenerates every figure from the committed code and cached sweep results, so
the paper's exhibits are reproducible artifacts rather than transient temp
files. Figures follow the numbering of the paper skeleton:

    F2  VGBF -- frequency/damping vs airspeed, QS and Peters, Fig-6 overlay
    F3  delta down-sweep branches (structural exhaustion)     [cached]
    F4  Barton verification bifurcation diagram               [cached]
    F5  NACA 0020 polar with the rig trim marked at Cl_max
    F6  cost-vs-constraint frontier, non-conservative region shaded
    F7  ONERA dynamic-stall subcritical branch                [cached]
    T1  results matrix (CSV + markdown)

F1 (schematic sub/supercritical) is drawn by hand. Cached inputs are the
JSONs written by stage0_barton_bifurcation.py, stage2b (delta sweep) and
stage4_onera.py; delete them to force a recompute there.

Run:  PYTHONPATH=src python examples/build_paper_figures.py
"""
import csv
import json
import os
import tempfile

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from mflco.model.michigan_params import (
    calibrate_michigan, section_from_params, structural_zeta)
from mflco.model.analysis import modal_analysis
from mflco.aero.quasi_steady import QuasiSteady
from mflco.aero.peters_finite import PetersFinite
from mflco.aero.qs_stall import load_polar

RIG_FLUTTER, RIG_FOLD = 13.19, 11.85
WINDOW_RIG = RIG_FLUTTER - RIG_FOLD               # 1.34 m/s
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "results")
TMP = tempfile.gettempdir()
os.makedirs(OUT, exist_ok=True)

cal = calibrate_michigan(zeta=structural_zeta())
sec = section_from_params(cal)


# --- F2: VGBF, both aero models ---------------------------------------------------
def f2_vgbf():
    Us = np.linspace(0.05, 5.0, 400)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
    for name, aero, col in (("QS", QuasiSteady(sec, 0.0), "tab:blue"),
                            ("Peters N=6", PetersFinite(sec, 0.0, N=6),
                             "tab:orange")):
        f, d, ums = [], [], []
        for u in Us:
            try:
                fr, dp = modal_analysis(sec, 0.0, u, aero)
            except ValueError:
                continue
            f.append(np.sort(fr)[:2])
            d.append(np.sort(dp)[:2])
            ums.append(cal.ustar_to_ms(u))
        f, d, ums = np.array(f), np.array(d), np.array(ums)
        hz = f * cal.omega_alpha / (2 * np.pi)
        for j in range(2):
            ax[0].plot(ums, hz[:, j], color=col, lw=1.4,
                       label=name if j == 0 else None)
            ax[1].plot(ums, d[:, j], color=col, lw=1.4,
                       label=name if j == 0 else None)
    # Garcia Perez Fig. 6 measured frequencies
    ax[0].plot([0, 4, 6, 8, 10, 11.8], [6.15, 6.13, 6.08, 6.02, 5.95, 5.90],
               "k^", ms=5, label="rig Fig.6 pitch")
    ax[0].plot([0, 4, 6, 8, 10, 11.8], [5.28, 5.33, 5.31, 5.34, 5.38, 5.45],
               "kv", ms=5, label="rig Fig.6 plunge")
    for a in ax:
        a.axvline(RIG_FLUTTER, color="r", ls="--", lw=1.1, label="rig flutter")
        a.axvline(RIG_FOLD, color="r", ls=":", lw=1.3, label="rig fold")
        a.set_xlabel("airspeed U [m/s]")
        a.grid(alpha=0.25)
        a.set_xlim(0, 16)
    ax[1].axhline(0.0, color="0.4", lw=0.9)
    ax[0].set_ylabel("frequency [Hz]"); ax[0].set_ylim(4.5, 7.0)
    ax[1].set_ylabel("damping ratio"); ax[1].set_ylim(-0.15, 0.25)
    ax[0].legend(fontsize=7.5); ax[1].legend(fontsize=7.5)
    fig.suptitle("F2  Velocity-frequency / velocity-damping: "
                 "wake lag recovers the QS flutter deficit", fontsize=11)
    plt.tight_layout()
    plt.savefig(f"{OUT}/F2_vgbf.png", dpi=170)
    plt.close()
    print("F2 -> results/F2_vgbf.png")


# --- F5: the polar, trim at Cl_max ------------------------------------------------
def f5_polar():
    ag, cl = load_polar()
    m = np.abs(ag) <= 30
    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    ax.plot(ag[m], cl[m], "o-", ms=3, color="tab:green",
            label="NACA 0020, Re 1.6e5 (Sheldahl-Klimas, interp.)")
    cl_trim = float(np.interp(10.0, ag, cl))
    ax.plot([10.0], [cl_trim], "r*", ms=16, zorder=5,
            label=f"rig trim 10 deg: Cl = {cl_trim:.3f} (at static Cl_max)")
    lin = np.linspace(-6, 12, 50)
    ax.plot(lin, 2 * np.pi * np.radians(lin), "k--", lw=1.0,
            label="attached-flow 2*pi*alpha")
    ax.annotate("local slope 0.63 /rad\n(90% collapse vs 2*pi)",
                xy=(10.0, cl_trim), xytext=(14.5, 0.35), fontsize=8.5,
                arrowprops=dict(arrowstyle="->", lw=0.9))
    ax.axvspan(10 - 14, 10 + 14, color="tab:red", alpha=0.07)
    ax.text(10, -0.55, "LCO incidence sweep at 14 deg amplitude",
            ha="center", fontsize=8, color="tab:red")
    ax.set_xlabel("angle of attack [deg]"); ax.set_ylabel("Cl")
    ax.set_title("F5  The rig trims at its airfoil's static Cl_max")
    ax.legend(fontsize=8, loc="upper left"); ax.grid(alpha=0.25)
    ax.set_xlim(-8, 30); ax.set_ylim(-0.7, 1.3)
    plt.tight_layout()
    plt.savefig(f"{OUT}/F5_polar_trim.png", dpi=170)
    plt.close()
    print("F5 -> results/F5_polar_trim.png")


# --- T1 + F6: results matrix and the frontier ------------------------------------
# constraint error is measured against the FOLD (11.85): a model that reports no
# fold hands the designer its Hopf, hence +10.95% into the unsafe region.
ROWS = [
    # model, class, U_flutter, fold?, U_fold_reported, cost/query s, executed
    ("QS", "inviscid-attached", 9.665, False, 9.665, 0.6, True),
    ("Peters N=6", "inviscid-attached", 13.148, False, 13.148, 1.3, True),
    ("QS+Stall", "viscous-static", None, False, None, 0.6, True),
    ("ONERA (Petot laws)", "viscous-dynamic", 10.875, True, 9.594, 28.3, True),
    ("UVLM", "inviscid-attached", None, None, None, None, False),
    ("Euler (steady + forced)", "inviscid ceiling", "--", False, "none",
     None, "ceiling executed; constraint query positioned"),
]


def t1_matrix():
    hdr = ["model", "class", "U_flutter [m/s]", "err vs 13.19 [%]", "fold?",
           "U_constraint [m/s]", "err vs 11.85 [%]", "cost/query [s]", "status"]
    lines = []
    for m, cls, uf, fold, ufold, cost, done in ROWS:
        e_f = (f"{100*(uf-RIG_FLUTTER)/RIG_FLUTTER:+.2f}"
               if isinstance(uf, float) else "--")
        e_c = (f"{100*(ufold-RIG_FOLD)/RIG_FOLD:+.2f}"
               if isinstance(ufold, float) else "--")
        uf_c = uf if isinstance(uf, str) else (f"{uf:.3f}" if uf else "none")
        ufold_c = ufold if isinstance(ufold, str) else (
            f"{ufold:.3f}" if ufold else "none")
        lines.append([m, cls, uf_c,
                      e_f, {True: "yes", False: "no", None: "--"}[fold],
                      ufold_c, e_c,
                      f"{cost:.1f}" if cost else "est.",
                      done if isinstance(done, str) else
                      ("executed" if done else "pending/positioned")])
    lines.append(["rig (experiment)", "--", f"{RIG_FLUTTER}", "--", "yes",
                  f"{RIG_FOLD}", "0.00", "--", "reference"])
    with open(f"{OUT}/T1_results_matrix.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(hdr); w.writerows(lines)
    with open(f"{OUT}/T1_results_matrix.md", "w") as f:
        f.write("| " + " | ".join(hdr) + " |\n")
        f.write("|" + "---|" * len(hdr) + "\n")
        for r in lines:
            f.write("| " + " | ".join(r) + " |\n")
    with open(f"{OUT}/T1_results_matrix.md", "a") as f:
        f.write("\n*QS reports its own flutter (9.665) as the boundary -- "
                "below the true fold, conservative only by accident of a "
                "-26.7% flutter error. QS+Stall predicts no instability at "
                "any airspeed: its reported boundary is unbounded -- "
                "qualitatively the least conservative entry if trusted. The "
                "+10.95% floor is the fold-blind class's best case, attained "
                "at converged flutter (Peters). Per-query costs: uniform "
                "settle primitive of Section 3.5, stamped on the evaluation "
                "machine.*\n")
    print("T1 -> results/T1_results_matrix.csv + .md")


def f6_frontier():
    fig, ax = plt.subplots(figsize=(8.4, 5.2))
    ax.axhspan(0, 20, color="tab:red", alpha=0.10)
    ax.axhline(0, color="k", lw=1.2)
    ax.text(60, 15.6, "NON-CONSERVATIVE: model reports a boundary\n"
            "ABOVE the true fold -- unsafe for design", fontsize=8.5, color="tab:red")
    ax.text(0.62, -4.2, "conservative (safe side)", fontsize=8.5,
            color="tab:green")
    floor = 100 * (RIG_FLUTTER - RIG_FOLD) / RIG_FOLD
    ax.axhline(floor, color="tab:red", ls="--", lw=1.3)
    ax.text(0.62, floor - 2.4, f"fold-blind class BEST CASE (converged "
            f"flutter): +{floor:.1f}%", fontsize=9, color="tab:red")
    ax.plot(1.3, floor, "o", ms=11, color="tab:orange", zorder=5)
    ax.annotate("Peters N=6\n(class best case)", (1.3, floor),
                textcoords="offset points", xytext=(10, 8), fontsize=8.5)
    ax.plot(0.6, -18.44, "o", ms=11, mfc="none", mec="tab:blue", mew=2,
            zorder=5)
    ax.annotate("QS: conservative BY ACCIDENT\n(flutter -26.7%: disqualified,"
                "\nnot informative)", (0.6, -18.44),
                textcoords="offset points", xytext=(-6, -34), fontsize=8)
    ax.plot(0.6, 18.6, "^", ms=11, mfc="none", mec="tab:green", mew=2,
            zorder=5)
    ax.annotate("QS+Stall: NO boundary reported\n(no instability predicted"
                " -> unbounded)", (0.6, 18.6), textcoords="offset points",
                xytext=(10, -8), fontsize=8, color="tab:green")
    ax.plot(28.3, -19.04, "o", ms=11, color="tab:purple", zorder=5)
    ax.annotate("ONERA (Petot laws): conservative\nBY MECHANISM -- right class,\nborrowed coefficients", (28.3, -19.04),
                textcoords="offset points", xytext=(16, 6),
                fontsize=8.5, color="tab:purple")
    for name, cost in (("Euler", 6710.0),):
        ax.plot(cost, floor, "s", ms=10, mfc="none", mec="0.5", ls="none")
        ax.annotate("Euler: 6,710 s per forced case (MEASURED);\n"
                    "constraint query ~ O(10) cases ~ 1e5 s\n"
                    "(fold-blind by constr.; ceiling executed, 5.3)",
                    (cost, floor),
                    textcoords="offset points", xytext=(-198, 24),
                    fontsize=8.5, color="0.4")
    ax.plot(0.0, 0.0, "r*", ms=18, zorder=6)
    ax.annotate("rig: true constraint\n(fold 11.85 m/s)", (0.0, 0.0),
                textcoords="offset points", xytext=(14, -26), fontsize=8.5)
    ax.set_xscale("symlog", linthresh=1.0)
    ax.set_xlabel("cost per fold-velocity query [s]  (log)")
    ax.set_ylabel("constraint error vs the fold [%]")
    ax.set_title("F6  Cost vs constraint accuracy: the fold-blind class scatters; "
                 "its best case is the floor")
    ax.grid(alpha=0.25); ax.set_ylim(-24, 20); ax.set_xlim(-0.3, 20000)
    plt.tight_layout()
    plt.savefig(f"{OUT}/F6_frontier.png", dpi=170)
    plt.close()
    print("F6 -> results/F6_frontier.png")


# --- F7: the ONERA subcritical branch ---------------------------------------------
def f7_onera():
    """ONERA branch: published Petot laws primary, r-law variant dashed."""
    nom = json.load(open(os.path.join(TMP, "stage4_onera.json")))["nominal"]
    fig, ax = plt.subplots(figsize=(8.4, 5.2))
    br = np.array(nom["branch"])
    ax.plot(br[:, 0], br[:, 1], "o-", ms=3.5, color="tab:purple",
            label="stable LCO branch (published NACA 0012 laws, TM-88917)")
    ax.plot([nom["U_fold"]], [nom["amp_at_fold"]], "v", color="tab:purple",
            ms=10, label=f"fold {nom['U_fold']:.2f} m/s "
                         f"(amp {nom['amp_at_fold']:.1f} deg)")
    hb = nom.get("hopf_bracket", [10.75, 11.0])
    ax.plot([0.5 * sum(hb)], [0.0], "*", color="tab:purple", ms=13,
            label=f"Hopf {0.5*sum(hb):.2f} m/s (bracket {hb[0]}-{hb[1]})")
    try:
        hyb = json.load(open(os.path.join(
            TMP, "stage4_onera_R2hybrid.json")))["nominal"]
        b2 = np.array(hyb["branch"])
        ax.plot(b2[:, 0], b2[:, 1], "s--", ms=2.5, color="tab:purple",
                alpha=0.45, label=f"r-law variant (R2=0.23): "
                                  f"fold {hyb['U_fold']:.2f}")
    except (FileNotFoundError, KeyError):
        pass
    ax.axvline(RIG_FLUTTER, color="r", ls="--", lw=1.1, label="rig flutter 13.19")
    ax.axvline(RIG_FOLD, color="r", ls=":", lw=1.4, label="rig fold 11.85")
    ax.axhspan(25, 45, color="0.5", alpha=0.08)
    ax.text(13.55, 40, "beyond static-polar validity\n(report with caution)",
            fontsize=8, color="0.35")
    ax.set_xlabel("airspeed U [m/s]")
    ax.set_ylabel("pitch LCO amplitude [deg]")
    ax.set_title("F7  ONERA dynamic stall: the subcritical branch the inviscid "
                 "axis could not produce")
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(alpha=0.25)
    ax.set_xlim(5.0, 16.5); ax.set_ylim(-1.5, 45)
    plt.tight_layout()
    plt.savefig(f"{OUT}/F7_onera_branch.png", dpi=170)
    plt.close()
    print("F7 -> results/F7_onera_branch.png")


# --- cached figures ----------------------------------------------------------------
def copy_cached():
    import shutil
    for src, dst in (("stage2b_delta_downsweep.png", "F3_delta_downsweep.png"),
                     ("stage0_barton_bifurcation.png",
                      "F4_barton_verification.png")):
        p = os.path.join(TMP, src)
        if os.path.exists(p):
            shutil.copy(p, f"{OUT}/{dst}")
            print(f"F{dst[1]} -> results/{dst}  (cached)")
        else:
            print(f"  [{dst} not cached -- run its stage script first]")


f2_vgbf()
f5_polar()
copy_cached()
t1_matrix()
f6_frontier()
f7_onera()
print(f"\nall artifacts -> {OUT}")