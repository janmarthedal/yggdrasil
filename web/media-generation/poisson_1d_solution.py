"""Generate poisson-1d-solution.svg.

Plots the analytical solution u(x) = x(1-x)/2 to the 1D Poisson problem
  -u'' = 1 on (0,1),  u(0) = u(1) = 0.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

out_path = Path(__file__).parent.parent / "media" / "poisson-1d-solution.svg"
out_path.parent.mkdir(parents=True, exist_ok=True)

x = np.linspace(0, 1, 300)
u = x * (1 - x) / 2

fig, ax = plt.subplots(figsize=(6, 3.5))

ax.plot(x, u, color="#1a5fa8", linewidth=2.5)
ax.fill_between(x, u, alpha=0.12, color="#1a5fa8")

# Boundary markers
ax.plot([0, 1], [0, 0], "o", color="#c0392b", markersize=7, zorder=3)

# Annotations
ax.annotate(r"$u(0)=0$", xy=(0, 0), xytext=(0.08, 0.035),
            fontsize=12, color="#c0392b",
            arrowprops=dict(arrowstyle="-|>", color="#c0392b", lw=1.2))
ax.annotate(r"$u(1)=0$", xy=(1, 0), xytext=(0.78, 0.035),
            fontsize=12, color="#c0392b",
            arrowprops=dict(arrowstyle="-|>", color="#c0392b", lw=1.2))
ax.text(0.5, u.max() + 0.005, r"$u(x)=\dfrac{x(1-x)}{2}$",
        fontsize=13, ha="center", va="bottom", color="#1a5fa8")

ax.set_xlabel(r"$x$", fontsize=13)
ax.set_ylabel(r"$u$", fontsize=13, rotation=0, labelpad=10)
ax.set_xlim(-0.05, 1.05)
ax.set_ylim(-0.01, 0.16)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()
fig.savefig(out_path, format="svg", bbox_inches="tight")
print(f"Saved {out_path}")
