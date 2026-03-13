"""Generate quad4-shape-functions.svg.

Shows the four bilinear shape functions of the Quad4 element as filled contour
plots over the reference square [0,1]^2.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

out_path = Path(__file__).parent.parent / "media" / "quad4-shape-functions.svg"
out_path.parent.mkdir(parents=True, exist_ok=True)

GRAY = "#555555"

n = 60
t = np.linspace(0, 1, n)
X, Y = np.meshgrid(t, t)

N = [
    (1.0 - X) * (1.0 - Y),
    X * (1.0 - Y),
    X * Y,
    (1.0 - X) * Y,
]
labels = [
    r"$N_0 = (1-\hat{x})(1-\hat{y})$",
    r"$N_1 = \hat{x}(1-\hat{y})$",
    r"$N_2 = \hat{x}\,\hat{y}$",
    r"$N_3 = (1-\hat{x})\,\hat{y}$",
]
node_x = [0.0, 1.0, 1.0, 0.0]
node_y = [0.0, 0.0, 1.0, 1.0]
node_labels = [r"$\hat{x}_0$", r"$\hat{x}_1$", r"$\hat{x}_2$", r"$\hat{x}_3$"]
offsets = [(-0.10, -0.10), (0.05, -0.10), (0.05, 0.05), (-0.10, 0.05)]

fig, axes = plt.subplots(1, 4, figsize=(11.5, 2.9))

for ax, Ni, label in zip(axes, N, labels):
    cf = ax.contourf(X, Y, Ni, levels=20, cmap="Blues")
    ax.contour(X, Y, Ni, levels=5, colors="#333333", linewidths=0.5, alpha=0.5)
    fig.colorbar(cf, ax=ax, fraction=0.05, pad=0.03)

    # Square outline
    sq = plt.Polygon([[0, 0], [1, 0], [1, 1], [0, 1]], closed=True,
                     edgecolor=GRAY, facecolor="none", linewidth=1.5, zorder=3)
    ax.add_patch(sq)

    # Node markers and labels
    ax.plot(node_x, node_y, "o", color=GRAY, markersize=6, zorder=5)
    for nx, ny, lbl, (ox, oy) in zip(node_x, node_y, node_labels, offsets):
        ax.text(nx + ox, ny + oy, lbl, fontsize=9, color=GRAY, ha="center")

    ax.set_xlim(-0.18, 1.13)
    ax.set_ylim(-0.18, 1.13)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(label, fontsize=9)

fig.tight_layout()
fig.savefig(out_path, format="svg", bbox_inches="tight")
print(f"Saved {out_path}")
