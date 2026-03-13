"""Generate quadrilateral-domain-quadrature.svg.

Shows tensor-product Gauss-Legendre quadrature points on [0,1]^2 for
order 7 (4x4 = 16 points). Marker size is proportional to the weight.
"""

from pathlib import Path

import matplotlib.pyplot as plt
from yggdrasil.refdomains.quadrilateral import QuadrilateralDomain

out_path = Path(__file__).parent.parent / "media" / "quadrilateral-domain-quadrature.svg"
out_path.parent.mkdir(parents=True, exist_ok=True)

BLUE = "#1a5fa8"

domain = QuadrilateralDomain()
pts, wts = domain.quadrature(7)
n = int(round(len(pts) ** 0.5))

fig, ax = plt.subplots(figsize=(3.5, 3.5))

ax.add_patch(plt.Polygon([[0, 0], [1, 0], [1, 1], [0, 1]], closed=True,
                          edgecolor="#555555", facecolor="#f0f4fa", linewidth=1.5))
ax.plot([0, 1, 1, 0], [0, 0, 1, 1], "o", color="#555555", markersize=5, zorder=3)

sizes = (wts / wts.max()) * 180
ax.scatter(pts[:, 0], pts[:, 1], s=sizes, color=BLUE, zorder=4)

ax.set_xlim(-0.08, 1.08)
ax.set_ylim(-0.08, 1.08)
ax.set_aspect("equal")
ax.axis("off")
ax.set_title(f"Order 7 ({n}×{n} = {n*n} points)", fontsize=11)

fig.tight_layout()
fig.savefig(out_path, format="svg", bbox_inches="tight")
print(f"Saved {out_path}")
