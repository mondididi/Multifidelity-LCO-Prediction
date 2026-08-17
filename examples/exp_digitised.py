"""Digitised experimental baselines for both wings -- single source of truth.

Rough manual digitisation (this work) from:

  BRISTOL:  Tartaruga et al. (2019, IFASD-155), Figure 6.
            Amplitude is HEAVE in mm (their diagram's DoF). Read tolerance
            ~ +/-0.15 m/s, +/-2 mm. Their text says the jump-down occurs
            "just above 16 m/s" and flutter "just after 24 m/s"; the plotted
            markers read ~17.15 and ~24.7 -- both readings are kept, the
            plotted ones are used.

  MICHIGAN: Garcia Perez et al. (2024), Fig. 9, PITCH panel (deg).
            Pitch chosen over plunge because the structural nonlinearity is
            the cubic PITCH spring, the model's amplitude variable and the
            beta anchor are pitch, so the overlay needs no unit conversion.
            Read tolerance ~ +/-0.05 m/s, +/-0.4 deg (source bars ~0.5 deg).

The experimental zero line is drawn as the subcritical picture:
    solid   0 -> fold        (equilibrium the only attractor)
    dotted  fold -> Hopf     (bistable: equilibrium coexists with the LCO)
    dashed  Hopf -> edge     (equilibrium unstable)

Import `draw_experiment(ax, wing)` for overlays, or run this file to write
results/exp_bristol_lco.csv, results/exp_michigan_lco.csv and a side-by-side
check figure against the sources (results/exp_digitised_check.png).

Run:  PYTHONPATH=src python examples/exp_digitised.py
"""
import csv

import numpy as np

BRISTOL = dict(
    units="heave_mm",
    fold=17.15,          # plotted jump marker (text: "just above 16")
    hopf=24.7,           # plotted star       (text: "just after 24")
    branch=[
        (17.15, 42), (19.0, 52), (20.0, 61), (20.7, 66), (21.4, 72),
        (22.6, 81), (24.0, 87), (25.0, 95), (26.0, 103), (27.0, 110)],
    unstable=[
        (17.3, 38), (18.0, 35), (19.0, 32), (20.0, 29), (21.0, 25),
        (22.0, 20), (23.0, 13), (24.0, 5), (24.7, 0)],
    eq=[14.5, 15.5, 16.5, 17.5, 18.0, 18.5, 19.0, 19.5, 20.0, 20.5, 21.0,
        21.5, 22.0, 22.5, 23.0, 23.5, 24.0, 24.5],
)

MICHIGAN = dict(
    units="pitch_deg",
    fold=11.85,
    hopf=13.19,
    branch=[
        (11.85, 3.9), (11.95, 4.15), (12.15, 4.4), (12.4, 5.7),
        (12.65, 9.6), (13.2, 14.1), (13.4, 14.6)],
    unstable=[],
    eq=[11.12, 11.43, 11.58, 12.15, 12.4, 12.65, 13.2],
)

WINGS = {"bristol": BRISTOL, "michigan": MICHIGAN}


def draw_experiment(ax, wing, x_max=None, color="k", label_prefix="exp",
                    show_branch=True):
    """Overlay one wing's digitised baseline in the subcritical-picture style."""
    d = WINGS[wing]
    x_max = x_max if x_max is not None else d["branch"][-1][0] + 0.6
    x_min = min(d["eq"]) - 0.4
    zlab = None if show_branch else f"{label_prefix} (digitised)"
    ax.plot([x_min, d["fold"]], [0, 0], "-", lw=1.6, color=color, label=zlab)
    ax.plot([d["fold"], d["hopf"]], [0, 0], ":", lw=1.6, color=color)
    ax.plot([d["hopf"], x_max], [0, 0], "--", lw=1.2, color=color, alpha=0.8)
    ax.plot(d["eq"], np.zeros(len(d["eq"])), ".", ms=6, color=color)
    if show_branch:
        u, a = zip(*d["branch"])
        ax.plot(u, a, "s-", ms=5, lw=1.3, mfc="none", mec=color, color=color,
                label=f"{label_prefix} (digitised)")
        if d["unstable"]:
            uu, ua = zip(*d["unstable"])
            ax.plot(uu, ua, "--", lw=1.2, color="0.45",
                    label=f"{label_prefix} sketched unstable br.")
    ax.plot([d["fold"]], [0], "x", ms=9, color=color)
    ax.plot([d["hopf"]], [0], "*", ms=12, color=color)
    return d


def write_csvs():
    for wing, d in WINGS.items():
        path = f"results/exp_{wing}_lco.csv"
        with open(path, "w", newline="") as f:
            f.write(f"# Digitised baseline ({d['units']}); provenance and "
                    "tolerance in examples/exp_digitised.py\n")
            w = csv.writer(f)
            w.writerow(["U_ms", "amplitude", "kind"])
            for U in d["eq"]:
                w.writerow([U, 0.0, "eq"])
            for U, a in d["branch"]:
                w.writerow([U, a, "lco"])
            for U, a in d["unstable"]:
                w.writerow([U, a, "unstable"])
            w.writerow([d["fold"], 0.0, "fold"])
            w.writerow([d["hopf"], 0.0, "hopf"])
        print(f"-> {path}")


if __name__ == "__main__":
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    write_csvs()
    fig, ax = plt.subplots(1, 2, figsize=(10.6, 3.9))
    draw_experiment(ax[0], "bristol", label_prefix="Bristol exp")
    ax[0].set_title("Bristol (heave, mm) -- check vs Tartaruga Fig. 6",
                    fontsize=9.5)
    ax[0].set_ylabel("heave amplitude [mm]")
    draw_experiment(ax[1], "michigan", label_prefix="Michigan exp")
    ax[1].set_title("Michigan (pitch, deg) -- check vs Garcia Perez Fig. 9",
                    fontsize=9.5)
    ax[1].set_ylabel("pitch amplitude [deg]")
    for a in ax:
        a.set_xlabel("airspeed U [m/s]")
        a.grid(alpha=0.25)
        a.legend(fontsize=7.5)
    plt.tight_layout()
    plt.savefig("results/exp_digitised_check.png", dpi=150)
    print("-> results/exp_digitised_check.png")