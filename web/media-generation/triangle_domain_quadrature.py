"""Generate triangle-domain-quadrature.svg.

Shows quadrature points on the reference triangle for orders 1 through 5
(1, 3, 4, 6, 7 points). Each panel depicts the reference triangle with
points scaled by |weight|. Positive weights are blue, negative are red.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

from yggdrasil.refdomains.triangle import TriangleDomain

out_path = Path(__file__).parent.parent / "media" / "triangle-domain-quadrature.svg"
out_path.parent.mkdir(parents=True, exist_ok=True)

BLUE = "#1a5fa8"
RED = "#c0392b"
GRAY = "#aaaaaa"

domain = TriangleDomain()
orders = [1, 2, 3, 4, 5]

# Reference triangle vertices
tri = plt.Polygon([[0, 0], [1, 0], [0, 1]], closed=True,
                  edgecolor="#555555", facecolor="#f0f4fa", linewidth=1.5)

fig, axes = plt.subplots(1, 5, figsize=(10, 2.6))

max_w = 0.5  # area of reference triangle — max possible weight

for ax, order in zip(axes, orders):
    pts, wts = domain.quadrature(order)
    n = len(pts)

    # Draw reference triangle
    ax.add_patch(plt.Polygon([[0, 0], [1, 0], [0, 1]], closed=True,
                             edgecolor="#555555", facecolor="#f0f4fa", linewidth=1.5))

    # Draw vertices
    ax.plot([0, 1, 0], [0, 0, 1], "o", color="#555555", markersize=4, zorder=3)

    # Draw quadrature points
    colors = [BLUE if w >= 0 else RED for w in wts]
    sizes = (np.abs(wts) / max_w) * 160
    ax.scatter(pts[:, 0], pts[:, 1], s=sizes, c=colors, zorder=4)

    ax.set_xlim(-0.08, 1.08)
    ax.set_ylim(-0.12, 1.08)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(f"order {order}\n({n} pt{'s' if n > 1 else ''})", fontsize=9)

# Legend
legend_handles = [
    mpatches.Patch(color=BLUE, label="positive weight"),
    mpatches.Patch(color=RED,  label="negative weight"),
]
fig.legend(handles=legend_handles, loc="lower center", ncol=2,
           fontsize=8, frameon=False, bbox_to_anchor=(0.5, -0.04))

fig.suptitle("Quadrature points on the reference triangle", fontsize=11, y=1.01)
fig.tight_layout()
fig.savefig(out_path, format="svg", bbox_inches="tight")
print(f"Saved {out_path}")
