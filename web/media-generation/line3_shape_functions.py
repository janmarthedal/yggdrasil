"""Generate line3-shape-functions.svg.

Shows the three quadratic Lagrange shape functions of the Line3 element
on [0, 1] with nodes at 0, 1, and 0.5.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

out_path = Path(__file__).parent.parent / "media" / "line3-shape-functions.svg"
out_path.parent.mkdir(parents=True, exist_ok=True)

BLUE = "#1a5fa8"
RED = "#c0392b"
GREEN = "#2a9d3e"
GRAY = "#555555"

x = np.linspace(0, 1, 300)
N0 = (1.0 - x) * (1.0 - 2.0 * x)
N1 = x * (2.0 * x - 1.0)
N2 = 4.0 * x * (1.0 - x)

fig, ax = plt.subplots(figsize=(5, 3.2))

ax.plot(x, N0, color=BLUE, linewidth=2,
        label=r"$N_0 = (1-\hat{x})(1-2\hat{x})$")
ax.plot(x, N1, color=RED, linewidth=2,
        label=r"$N_1 = \hat{x}(2\hat{x}-1)$")
ax.plot(x, N2, color=GREEN, linewidth=2,
        label=r"$N_2 = 4\hat{x}(1-\hat{x})$")

node_x = [0.0, 1.0, 0.5]
ax.plot(node_x, [0, 0, 0], "|", color=GRAY, markersize=14,
        markeredgewidth=2, zorder=5)
ax.plot(node_x, [0, 0, 0], "o", color=GRAY, markersize=6, zorder=5)

ax.axhline(0, color="#aaaaaa", linewidth=0.8, zorder=0)
ax.set_xlim(-0.05, 1.05)
ax.set_ylim(-0.35, 1.15)
ax.set_xlabel(r"$\hat{x}$")
ax.legend(fontsize=10)
ax.set_title("Line3 shape functions", fontsize=12)

fig.tight_layout()
fig.savefig(out_path, format="svg", bbox_inches="tight")
print(f"Saved {out_path}")
