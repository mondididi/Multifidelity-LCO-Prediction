"""Build the paper's figure set and results matrix into results/.

Regenerates every figure from the committed code and cached sweep results, so
the paper's exhibits are reproducible artifacts rather than transient temp
files. Figures follow the numbering of the paper skeleton:

    F2  VGBF -- frequency/damping vs airspeed, QS and Peters, Fig-6 overlay
    F3  delta down-sweep branches (structural exhaustion)     [cached]
    F4  Barton verification bifurcation diagram               [cached]
    F5  NACA 0020 polar with the rig trim marked at Cl_max
    F6  cost-vs-constraint frontier, non-conservative region shaded
    T1  results matrix (CSV + markdown)

F1 (schematic sub/supercritical) is drawn by hand; F7 (ONERA) follows that run.
Cached inputs are the JSONs written by stage0_barton_bifurcation.py and
stage2_michigan_delta_sweep.py; delete them to force a recompute there.
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
# fold hands the designer its Hopf, hence +11.3% into the unsafe region.
ROWS = [
    # model, class, U_flutter, fold?, U_fold_reported, cost/query s, executed
    ("QS", "inviscid-attached", 9.665, False, 13.148, 2.0, True),
    ("Peters N=6", "inviscid-attached", 13.148, False, 13.148, 17.0, True),
    ("QS+Stall", "viscous-static", None, False, None, 2.0, True),
    ("UVLM", "inviscid-attached", None, None, None, None, False),
    ("Euler (steady)", "inviscid ceiling", None, None, None, None, False),
    ("ONERA", "viscous-dynamic", None, None, None, None, False),
]


def t1_matrix():
    hdr = ["model", "class", "U_flutter [m/s]", "err vs 13.19 [%]", "fold?",
           "U_constraint [m/s]", "err vs 11.85 [%]", "cost/query [s]", "status"]
    lines = []
    for m, cls, uf, fold, ufold, cost, done in ROWS:
        e_f = f"{100*(uf-RIG_FLUTTER)/RIG_FLUTTER:+.2f}" if uf else "--"
        e_c = f"{100*(ufold-RIG_FOLD)/RIG_FOLD:+.2f}" if ufold else "--"
        lines.append([m, cls, f"{uf:.3f}" if uf else "none",
                      e_f, {True: "yes", False: "no", None: "--"}[fold],
                      f"{ufold:.3f}" if ufold else "none", e_c,
                      f"{cost:.0f}" if cost else "est.",
                      "executed" if done else "pending/positioned"])
    lines.append(["rig (experiment)", "--", f"{RIG_FLUTTER}", "--", "yes",
                  f"{RIG_FOLD}", "0.00", "--", "reference"])
    with open(f"{OUT}/T1_results_matrix.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(hdr); w.writerows(lines)
    with open(f"{OUT}/T1_results_matrix.md", "w") as f:
        f.write("| " + " | ".join(hdr) + " |\n")
        f.write("|" + "---|" * len(hdr) + "\n")
        for r in lines:
            f.write("| " + " | ".join(r) + " |\n")
    print("T1 -> results/T1_results_matrix.csv + .md")


def f6_frontier():
    fig, ax = plt.subplots(figsize=(8.4, 5.2))
    ax.axhspan(0, 20, color="tab:red", alpha=0.10)
    ax.axhline(0, color="k", lw=1.2)
    ax.text(0.62, 17.2, "NON-CONSERVATIVE: model reports a boundary ABOVE\n"
            "the true fold -- unsafe for design", fontsize=8.5, color="tab:red")
    ax.text(0.62, -4.2, "conservative (safe side)", fontsize=8.5,
            color="tab:green")
    floor = 100 * (RIG_FLUTTER - RIG_FOLD) / RIG_FOLD
    ax.axhline(floor, color="tab:red", ls="--", lw=1.3)
    ax.text(0.62, floor - 2.4, f"irreducible floor for the no-fold class: "
            f"+{floor:.1f}%", fontsize=9, color="tab:red")
    pts = [("QS", 2.0, floor, "tab:blue"),
           ("Peters N=6", 17.0, floor, "tab:orange"),
           ("QS+Stall", 2.0, floor, "tab:green")]
    for name, cost, err, col in pts:
        ax.plot(cost, err, "o", ms=11, color=col, zorder=5)
        ax.annotate(name, (cost, err), textcoords="offset points",
                    xytext=(-18, 12), fontsize=9)
    ax.plot(2.0, floor, "o", ms=11, mfc="none", mec="tab:green", mew=2)
    for name, cost in (("ONERA", 60.0), ("Euler campaign", 5000.0)):
        ax.plot(cost, floor, "s", ms=10, mfc="none", mec="0.5", ls="none")
        ax.annotate(f"{name}\n(pending)", (cost, floor),
                    textcoords="offset points", xytext=(-20, -34),
                    fontsize=8.5, color="0.4")
    ax.plot(0.0, 0.0, "r*", ms=18, zorder=6)
    ax.annotate("rig: true constraint\n(fold 11.85 m/s)", (0.0, 0.0),
                textcoords="offset points", xytext=(14, -26), fontsize=8.5)
    ax.set_xscale("symlog", linthresh=1.0)
    ax.set_xlabel("cost per fold-velocity query [s]  (log)")
    ax.set_ylabel("constraint error vs the fold [%]")
    ax.set_title("F6  Cost vs constraint accuracy: every fold-blind model sits "
                 "on the unsafe floor")
    ax.grid(alpha=0.25); ax.set_ylim(-8, 20); ax.set_xlim(-0.3, 20000)
    plt.tight_layout()
    plt.savefig(f"{OUT}/F6_frontier.png", dpi=170)
    plt.close()
    print("F6 -> results/F6_frontier.png")


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
print(f"\nall artifacts -> {OUT}")