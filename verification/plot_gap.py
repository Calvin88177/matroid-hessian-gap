import csv, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from collections import defaultdict

rows = list(csv.DictReader(open("gap_normalized.csv")))
order = ["complete K_n", "bipartite K_{m,m}", "hypercube Q_d", "square grid m x m", "triangular grid (m,m)", "cycle C_n"]
label = {"complete K_n": "Kₙ", "bipartite K_{m,m}": "Kₘ,ₘ", "hypercube Q_d": "Hypercube",
         "square grid m x m": "Square grid", "triangular grid (m,m)": "Triangular grid", "cycle C_n": "Cycle"}
colors = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
markers = ["o", "s", "D", "^", "v", "P"]
data = defaultdict(list)
for r in rows:
    if int(r["edges"]) >= 4:
        data[r["family"]].append((int(r["edges"]), float(r["gamma_uniform"]), float(r["gamma_min_ratio2"])))

ink, muted, grid, surf = "#0b0b0b", "#898781", "#e1e0d9", "#fcfcfb"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})
fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), sharey=True, facecolor=surf)
panels = [(1, "Uniform weights (x = 1)"), (2, "Weights within a factor of 2 (worst of 20 samples)")]
for ax, (col, title) in zip(axes, panels):
    ax.set_facecolor(surf)
    for fam, c, m in zip(order, colors, markers):
        pts = sorted(data[fam])
        xs = [p[0] for p in pts]; ys = [p[col] for p in pts]
        ax.plot(xs, ys, color=c, lw=2, marker=m, ms=5, markeredgecolor=surf, markeredgewidth=1, label=label[fam])
    ax.set_xscale("log"); ax.set_ylim(0, 1.05)
    ax.set_title(title, loc="left", color=ink, fontsize=11)
    ax.set_xlabel("Number of edges |E| (log scale)", color=muted)
    ax.grid(True, axis="y", color=grid, lw=0.8); ax.grid(False, axis="x")
    for s in ["top", "right"]: ax.spines[s].set_visible(False)
    for s in ["left", "bottom"]: ax.spines[s].set_color("#c3c2b7")
    ax.tick_params(colors=muted)
axes[0].set_ylabel("Normalized gap  γ = (n−1)·min|λ₋| / λ₊", color=ink)
axes[1].legend(frameon=False, loc="upper right", fontsize=9, labelcolor=ink, ncol=2)
fig.suptitle("The Hessian's normalized gap does not decay with graph size", x=0.01, ha="left", color=ink, fontsize=13, fontweight="bold")
fig.tight_layout(rect=(0, 0, 1, 0.94))
fig.savefig("../figures/gap_vs_size.png", dpi=160, facecolor=surf)
print("saved")
