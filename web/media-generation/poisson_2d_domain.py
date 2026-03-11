"""Generate poisson-2d-domain.svg.

Illustrates the domain Ω = [-1, 1]² minus a circular hole,
with Γ_D (square sides, Dirichlet) and Γ_N (circle, Neumann) labelled.
"""

from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.path import Path as MPath

# Domain parameters (matching poisson_hole.py)
HALF = 1.0
CX, CY, R = 0.3, 0.2, 0.5

out_path = Path(__file__).parent.parent / "media" / "poisson-2d-domain.svg"
out_path.parent.mkdir(parents=True, exist_ok=True)

fig, ax = plt.subplots(figsize=(5, 5))

# --- Domain fill (square with circular hole) using a compound path ---
# Square vertices (counter-clockwise)
sq = np.array([[-HALF, -HALF], [HALF, -HALF], [HALF, HALF], [-HALF, HALF], [-HALF, -HALF]])
# Circle vertices (clockwise, so the interior is excluded)
theta = np.linspace(0, 2 * np.pi, 128)
circle = np.column_stack([CX + R * np.cos(theta), CY + R * np.sin(theta)])[::-1]
circle = np.vstack([circle, circle[0]])

verts = np.vstack([sq, circle])
codes = (
    [MPath.MOVETO] + [MPath.LINETO] * 3 + [MPath.CLOSEPOLY]
    + [MPath.MOVETO] + [MPath.LINETO] * (len(circle) - 2) + [MPath.CLOSEPOLY]
)
domain_path = MPath(verts, codes)
patch = mpatches.PathPatch(domain_path, facecolor="#dce8f5", edgecolor="none", zorder=1)
ax.add_patch(patch)

# --- Boundaries ---
# Γ_D: square sides (blue)
square_x = [-HALF, HALF, HALF, -HALF, -HALF]
square_y = [-HALF, -HALF, HALF, HALF, -HALF]
ax.plot(square_x, square_y, color="#1a5fa8", linewidth=2.5, zorder=2)

# Γ_N: circular boundary (red)
theta_plot = np.linspace(0, 2 * np.pi, 256)
ax.plot(CX + R * np.cos(theta_plot), CY + R * np.sin(theta_plot),
        color="#c0392b", linewidth=2.5, zorder=2)

# --- Labels ---
# Ω in the interior (away from the hole)
ax.text(-0.65, -0.65, r"$\Omega$", fontsize=22, ha="center", va="center",
        color="#222", zorder=3)

# Γ_D on the top side
ax.annotate(r"$\Gamma_D$", xy=(0.5, HALF), xytext=(0.75, 1.35),
            fontsize=16, color="#1a5fa8",
            arrowprops=dict(arrowstyle="-|>", color="#1a5fa8", lw=1.4),
            ha="center", zorder=3)

# Γ_N on the circle (point at top of circle)
circle_label_pt = (CX + R * np.cos(np.pi * 0.65), CY + R * np.sin(np.pi * 0.65))
ax.annotate(r"$\Gamma_N$", xy=circle_label_pt, xytext=(-0.6, 0.75),
            fontsize=16, color="#c0392b",
            arrowprops=dict(arrowstyle="-|>", color="#c0392b", lw=1.4),
            ha="center", zorder=3)

# --- Axes ---
ax.set_xlim(-1.5, 1.5)
ax.set_ylim(-1.5, 1.5)
ax.set_aspect("equal")
ax.axis("off")

plt.tight_layout()
fig.savefig(out_path, format="svg", bbox_inches="tight")
print(f"Saved {out_path}")
