"""Figure: the cost-constraint frontier (Michigan wing, four rungs).

Each rung is placed by the machine-stamped cost of one settled query
(examples/time_costs.py) and by the error of THE BOUNDARY IT ACTUALLY
REPORTS -- its fold if it has one, otherwise its own flutter speed --
against the measured fold at 11.85 m/s (values from the stage-3/4/5
records). The fold-blind class scatters: QS lands conservative by
accident (-18.4%), Peters attains the class best case (+10.95%, one
refinement short of the +11.3% perfect-flutter floor), the static-stall
rung reports NO boundary at all. Only the ONERA rung is conservative by
mechanism (-19.0%). Euler appears as future work: cost only, no boundary.

Writes: results/F6_frontier.png
Run:    PYTHONPATH=src python examples/fig_frontier.py   (seconds)
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# (cost s/query, reported-boundary error %, label, colour)   -- provenance:
# costs = examples/time_costs.py; boundaries = stage-3/4/5 records
DATA = [
    (0.6, -18.4, "QS\n(own flutter; accidental)", "tab:blue"),
    (1.3, +10.95, "Peters N=6\n(own flutter; class best)", "tab:orange"),
    (28.3, -19.0, "ONERA\n(fold; conservative by mechanism)", "tab:purple"),
]
UNBOUNDED = (0.6, "QS+Stall: no boundary reported")
EULER = (6710.0, "Euler (future work):\n6.7e3 s/case, no boundary")

fig, ax = plt.subplots(figsize=(7.4, 4.5))
ax.axhline(0, color="0.35", lw=1.0)
ax.axhline(11.3, color="tab:red", ls="--", lw=1.3)
ax.text(0.42, 11.9, "perfect-flutter floor of the fold-blind class (+11.3%)",
        color="tab:red", fontsize=8)
for cost, err, lab, c in DATA:
    ax.plot([cost], [err], "o", ms=10, color=c)
    ax.annotate(lab, (cost, err), textcoords="offset points", xytext=(10, -4),
                fontsize=8, color=c)
ax.plot([UNBOUNDED[0]], [24.5], "^", ms=11, mfc="none", mec="tab:green",
        mew=2)
ax.annotate(UNBOUNDED[1] + "\n(unbounded)", (UNBOUNDED[0], 24.5),
            textcoords="offset points", xytext=(10, -6), fontsize=8,
            color="tab:green")
ax.plot([EULER[0]], [0], "s", ms=10, mfc="none", mec="0.5", mew=1.6)
ax.annotate(EULER[1], (EULER[0], 0), textcoords="offset points",
            xytext=(-92, 12), fontsize=7.5, color="0.4")
ax.set_xscale("log")
ax.set_xlim(0.3, 2e4)
ax.set_ylim(-27, 29)
ax.set_xlabel("cost of one settled boundary query [s]  (log)")
ax.set_ylabel("reported-boundary error vs measured fold [%]")
ax.set_title("Michigan wing: cost vs constraint accuracy -- flutter accuracy "
             "buys nothing inside the fold-blind class", fontsize=9.5)
ax.grid(alpha=0.25, which="both")
plt.tight_layout()
plt.savefig("results/F6_frontier.png", dpi=170)
print("-> results/F6_frontier.png")