"""Generate line2-shape-functions.svg.

Shows the two linear Lagrange shape functions N0 = 1-x and N1 = x
of the Line2 element on [0, 1].
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

out_path = Path(__file__).parent.parent / "media" / "line2-shape-functions.svg"
out_path.parent.mkdir(parents=True, exist_ok=True)

BLUE = "#1a5fa8"
RED = "#c0392b"
GRAY = "#555555"

x = np.linspace(0, 1, 300)
N0 = 1.0 - x
N1 = x

fig, ax = plt.subplots(figsize=(5, 3))

ax.plot(x, N0, color=BLUE, linewidth=2, label=r"$N_0 = 1 - \hat{x}$")
ax.plot(x, N1, color=RED, linewidth=2, label=r"$N_1 = \hat{x}$")

ax.plot([0, 1], [0, 0], "|", color=GRAY, markersize=14, markeredgewidth=2, zorder=5)
ax.plot([0, 1], [0, 0], "o", color=GRAY, markersize=6, zorder=5)

ax.axhline(0, color="#aaaaaa", linewidth=0.8, zorder=0)
ax.set_xlim(-0.05, 1.05)
ax.set_ylim(-0.15, 1.15)
ax.set_xlabel(r"$\hat{x}$")
ax.legend(fontsize=11)
ax.set_title("Line2 shape functions", fontsize=12)

fig.tight_layout()
fig.savefig(out_path, format="svg", bbox_inches="tight")
print(f"Saved {out_path}")
