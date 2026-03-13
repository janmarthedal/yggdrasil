"""Generate line-domain-quadrature.svg.

Shows Gauss-Legendre quadrature points on [0, 1] for orders 1, 3, 5, 7
(i.e. 1, 2, 3, 4 quadrature points). Each row depicts the reference
interval with points scaled by weight.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from yggdrasil.refdomains.line import LineDomain

out_path = Path(__file__).parent.parent / "media" / "line-domain-quadrature.svg"
out_path.parent.mkdir(parents=True, exist_ok=True)

BLUE = "#1a5fa8"

domain = LineDomain()

# Orders to display (polynomial degree integrated exactly)
orders = [1, 3, 5, 7]

rows = len(orders)
fig, axes = plt.subplots(rows, 1, figsize=(7, 3.2))

for ax, order in zip(axes, orders):
    pts_arr, wts = domain.quadrature(order)
    pts = pts_arr[:, 0]
    n = len(pts)

    # Draw the reference interval
    ax.hlines(0, 0, 1, colors="#aaaaaa", linewidth=2, zorder=1)
    ax.plot([0, 1], [0, 0], "|", color="#555555", markersize=12, markeredgewidth=2)

    # Draw quadrature points; size encodes weight
    max_w = 0.5  # max weight on [0,1] is 0.5 (n=1)
    sizes = (wts / max_w) * 180
    ax.scatter(pts, np.zeros_like(pts), s=sizes, color=BLUE, zorder=3)

    # Annotate each point with its coordinate
    for x, w in zip(pts, wts):
        ax.text(x, 0.18, f"{x:.4f}", ha="center", va="bottom",
                fontsize=8, color=BLUE)

    # Row label on the left
    ax.text(-0.03, 0, f"$n={n}$", ha="right", va="center",
            fontsize=11, color="#333333")

    ax.set_xlim(-0.12, 1.12)
    ax.set_ylim(-0.5, 0.55)
    ax.axis("off")

axes[0].set_title("Gauss–Legendre points on $[0, 1]$", fontsize=12, pad=8)

fig.tight_layout(h_pad=0.2)
fig.savefig(out_path, format="svg", bbox_inches="tight")
print(f"Saved {out_path}")
