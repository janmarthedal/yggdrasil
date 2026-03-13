"""Generate tri3-shape-functions.svg.

Shows the three linear shape functions (barycentric coordinates) of the
Tri3 element as filled contour plots over the reference triangle.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np

out_path = Path(__file__).parent.parent / "media" / "tri3-shape-functions.svg"
out_path.parent.mkdir(parents=True, exist_ok=True)

GRAY = "#555555"

# Sample points inside the reference triangle
n = 60
xi_vals = np.linspace(0, 1, n)
pts = [(x, y) for x in xi_vals for y in xi_vals if x + y <= 1.0]
pts = np.array(pts)
x, y = pts[:, 0], pts[:, 1]

N = [1.0 - x - y, x, y]
labels = [r"$N_0 = 1 - \hat{x} - \hat{y}$", r"$N_1 = \hat{x}$", r"$N_2 = \hat{y}$"]
node_x = [0.0, 1.0, 0.0]
node_y = [0.0, 0.0, 1.0]
node_labels = [r"$\hat{x}_0$", r"$\hat{x}_1$", r"$\hat{x}_2$"]
offsets = [(-0.09, -0.09), (0.04, -0.09), (-0.09, 0.04)]

triang = mtri.Triangulation(x, y)

fig, axes = plt.subplots(1, 3, figsize=(8.5, 2.9))

for ax, Ni, label in zip(axes, N, labels):
    tcf = ax.tricontourf(triang, Ni, levels=20, cmap="Blues")
    ax.tricontour(triang, Ni, levels=5, colors="#333333", linewidths=0.5, alpha=0.5)
    fig.colorbar(tcf, ax=ax, fraction=0.05, pad=0.03)

    # Triangle outline
    tri_patch = plt.Polygon([[0, 0], [1, 0], [0, 1]], closed=True,
                             edgecolor=GRAY, facecolor="none", linewidth=1.5, zorder=3)
    ax.add_patch(tri_patch)

    # Node markers and labels
    ax.plot(node_x, node_y, "o", color=GRAY, markersize=6, zorder=5)
    for nx, ny, lbl, (ox, oy) in zip(node_x, node_y, node_labels, offsets):
        ax.text(nx + ox, ny + oy, lbl, fontsize=9, color=GRAY, ha="center")

    ax.set_xlim(-0.15, 1.1)
    ax.set_ylim(-0.15, 1.1)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(label, fontsize=10)

fig.tight_layout()
fig.savefig(out_path, format="svg", bbox_inches="tight")
print(f"Saved {out_path}")
