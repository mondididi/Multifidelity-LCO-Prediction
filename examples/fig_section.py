"""Figure: the 2-DOF pitch-plunge typical section (the Section 2.1 schematic).

Pure schematic -- no data dependencies. Names every point Sections 2.1-2.2
refer to: the elastic axis (where all structural attachment acts), the centre
of mass at x_alpha*b aft of it, mid-chord (the datum for a), the quarter-chord
(where the circulatory lift acts) and the three-quarter-chord (where alpha_eff
is sampled). The only nonlinearity in the model sits in the pitch spring, and
the figure says so.

Drawn generically: the elastic axis is forward of mid-chord (a < 0, as on both
wings) but is NOT drawn on top of the quarter-chord, so the two labels stay
legible. Both wings actually have a = -1/2 exactly -- that belongs in the
caption, not baked into the drawing.

Writes: results/F2_section.png and results/F2_section.pdf
Run:    PYTHONPATH=src python examples/fig_section.py   (seconds)
"""
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Polygon

# ---- geometry, in semi-chords with mid-chord at the origin -----------------------
LE, TE = -1.0, 1.0          # leading / trailing edge
C4, C34 = -0.5, 0.5         # quarter-chord, three-quarter-chord
A_EA = -0.30                # elastic axis as drawn (a < 0: forward of mid-chord)
X_CG = 0.16                 # centre of mass, x_alpha*b aft of the EA
BRACKET = -0.36             # underside of the mounting bracket
GROUND = -1.58              # support datum
THICK = 0.18                # NACA 0018-ish profile; both rigs are symmetric

INK = "0.15"
ACC = "#1f6fb4"             # structural attachment
AERO = "#c0392b"            # aerodynamic reference points and loads
DIM = "0.45"


def lead(color="0.55"):
    return dict(arrowstyle="-", color=color, lw=0.9, shrinkA=1, shrinkB=3)


def naca_symmetric(t=THICK, n=200):
    """Half-thickness of a symmetric 4-digit section, chord 2 (LE at -1)."""
    xc = np.linspace(0.0, 1.0, n)
    yt = 5 * t * (0.2969 * np.sqrt(xc) - 0.1260 * xc - 0.3516 * xc ** 2
                  + 0.2843 * xc ** 3 - 0.1015 * xc ** 4)
    return 2 * xc - 1.0, 2 * yt


def coil(ax, x, y0, y1, n=7, w=0.085, **kw):
    ys = np.linspace(y0, y1, 2 * n + 3)
    xs = np.full_like(ys, x)
    xs[1:-1:2] = x - w
    xs[2:-1:2] = x + w
    ax.plot(xs, ys, lw=1.5, solid_joinstyle="miter", **kw)


def dashpot(ax, x, y0, y1, w=0.085, **kw):
    ym = 0.5 * (y0 + y1)
    ax.plot([x, x], [y0, ym + 0.13], lw=1.5, **kw)
    ax.plot([x - w, x - w, x + w, x + w],
            [ym + 0.13, ym - 0.09, ym - 0.09, ym + 0.13], lw=1.5, **kw)
    ax.plot([x - w * 0.78, x + w * 0.78], [ym + 0.02, ym + 0.02], lw=2.8, **kw)
    ax.plot([x, x], [ym + 0.02, y1], lw=1.5, **kw)


def spiral(ax, x, y, turns=2.2, r0=0.035, r1=0.135, n=240, **kw):
    th = np.linspace(0.0, 2 * np.pi * turns, n)
    r = np.linspace(r0, r1, n)
    ax.plot(x + r * np.cos(th), y + r * np.sin(th), lw=1.4, **kw)


def hatched_ground(ax, x0, x1, y, n=13, h=0.11):
    ax.plot([x0, x1], [y, y], color=INK, lw=1.8)
    for xs in np.linspace(x0, x1, n)[:-1]:
        ax.plot([xs, xs + h], [y, y - h], color=INK, lw=1.0)


def dim(ax, x0, x1, y, label, dy=0.15, tick_to=None):
    ax.annotate("", xy=(x0, y), xytext=(x1, y),
                arrowprops=dict(arrowstyle="<->", color=DIM, lw=1.0))
    ax.text(0.5 * (x0 + x1), y + dy, label, fontsize=10, color="0.25",
            ha="center", va="center")
    if tick_to is not None:
        for xv in (x0, x1):
            ax.plot([xv, xv], [y + 0.05, tick_to], color="0.75", lw=0.7,
                    zorder=1)


fig, ax = plt.subplots(figsize=(9.0, 6.6))

# ---- the section -----------------------------------------------------------------
xf, yf = naca_symmetric()
ax.add_patch(Polygon(np.column_stack([np.r_[xf, xf[::-1]], np.r_[yf, -yf[::-1]]]),
                     closed=True, facecolor="0.92", edgecolor=INK, lw=1.7,
                     zorder=3))
ax.plot([LE, TE], [0, 0], color="0.45", lw=0.9, ls=(0, (6, 4)), zorder=4)

# ---- mounting: plunge spring + damper, torsional spring at the elastic axis -------
ax.plot([A_EA, A_EA], [-0.12, BRACKET], color=ACC, lw=1.6, zorder=2)
ax.plot([A_EA - 0.32, A_EA + 0.32], [BRACKET, BRACKET], color=ACC, lw=1.6)
coil(ax, A_EA - 0.32, BRACKET, GROUND, color=ACC)
dashpot(ax, A_EA + 0.32, BRACKET, GROUND, color=ACC)
hatched_ground(ax, A_EA - 0.64, A_EA + 0.64, GROUND)
spiral(ax, A_EA, 0.0, color=ACC, zorder=5)

ax.text(A_EA - 0.46, 0.5 * (BRACKET + GROUND), r"$k_\xi$", color=ACC,
        fontsize=12, ha="right", va="center")
ax.text(A_EA + 0.46, 0.5 * (BRACKET + GROUND), r"$c_\xi$", color=ACC,
        fontsize=12, ha="left", va="center")

# ---- named points ----------------------------------------------------------------
ax.plot([C4], [0], "s", color=AERO, ms=7, zorder=7)            # quarter-chord
ax.plot([C34], [0], "^", color=AERO, ms=8, zorder=7)           # three-quarter
ax.plot([0.0], [0], "|", color=INK, ms=11, mew=1.5, zorder=7)  # mid-chord
ax.plot([A_EA], [0], "o", color=ACC, ms=9, zorder=8)           # elastic axis
ax.plot([X_CG], [0], "o", color=INK, ms=10, mfc="white", mew=1.6, zorder=8)
ax.plot([X_CG], [0], "o", color=INK, ms=3.4, zorder=9)         # centre of mass

ax.annotate("elastic axis\n(all structural attachment)",
            xy=(A_EA - 0.04, -0.06), xytext=(-2.70, -0.66), fontsize=10.5,
            color=ACC, ha="left", va="center", arrowprops=lead(ACC))
ax.annotate(r"pitch spring $k_\alpha(\alpha)$: the model's" "\n"
            r"only nonlinearity; damper $c_\alpha$",
            xy=(A_EA + 0.15, 0.09), xytext=(1.22, 1.06), fontsize=10.5,
            color=ACC, ha="left", va="center", arrowprops=lead(ACC))
ax.annotate("centre of mass", xy=(X_CG + 0.04, -0.06), xytext=(1.22, -0.62),
            fontsize=10.5, color=INK, ha="left", va="center",
            arrowprops=lead())
ax.annotate(r"three-quarter-chord:" "\n" r"$\alpha_{\rm eff}$ sampled here",
            xy=(C34 + 0.03, 0.09), xytext=(1.22, 0.58), fontsize=10.5,
            color=AERO, ha="left", va="center", arrowprops=lead(AERO))
ax.text(0.0, 0.32, "mid-chord", fontsize=10, color="0.30", ha="center")

# ---- degrees of freedom ----------------------------------------------------------
ax.add_patch(FancyArrowPatch((-1.58, 0.24), (-1.58, -0.36), arrowstyle="-|>",
                             mutation_scale=14, color=INK, lw=1.7))
ax.text(-1.68, -0.02, r"$h$", fontsize=13, ha="right", va="center")
ax.text(-1.68, -0.26, "(+ down)", fontsize=9, ha="right", va="center",
        color="0.35")

ax.add_patch(FancyArrowPatch((-1.10, 0.02), (-0.74, 0.44), arrowstyle="-|>",
                             mutation_scale=14, color=INK, lw=1.7,
                             connectionstyle="arc3,rad=0.38"))
ax.text(-1.10, 0.54, r"$\alpha$ (+ nose-up)", fontsize=11.5, ha="right")

# ---- freestream and lift ---------------------------------------------------------
for yy in (0.16, -0.16):
    ax.add_patch(FancyArrowPatch((-2.70, yy), (-2.14, yy), arrowstyle="-|>",
                                 mutation_scale=12, color="0.5", lw=1.4))
ax.text(-2.42, 0.34, r"$U$", fontsize=13, ha="center", color="0.35")

ax.add_patch(FancyArrowPatch((C4, 0.12), (C4, 0.54), arrowstyle="-|>",
                             mutation_scale=14, color=AERO, lw=1.8))
ax.text(C4 - 0.02, 0.64, r"$L$ (quarter-chord)", fontsize=11, color=AERO,
        ha="center")

# ---- dimensions ------------------------------------------------------------------
dim(ax, LE, 0.0, 1.66, r"$b$", tick_to=0.50)
dim(ax, 0.0, TE, 1.66, r"$b$", tick_to=0.50)
dim(ax, A_EA, X_CG, 1.28, r"$x_\alpha b$", tick_to=0.50)
dim(ax, A_EA, 0.0, 0.92, r"$ab$")

ax.set_xlim(-2.78, 2.72)
ax.set_ylim(-1.92, 1.92)
ax.set_aspect("equal")
ax.axis("off")
plt.savefig("results/F2_section.png", dpi=230, bbox_inches="tight")
plt.savefig("results/F2_section.pdf", bbox_inches="tight")
print("-> results/F2_section.png / .pdf")