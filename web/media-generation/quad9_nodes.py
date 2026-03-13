"""Generate quad9-nodes.svg.

Shows the reference square with all 9 nodes of the Quad9 element
labelled by index and coordinate.
"""

from pathlib import Path

import matplotlib.pyplot as plt

out_path = Path(__file__).parent.parent / "media" / "quad9-nodes.svg"
out_path.parent.mkdir(parents=True, exist_ok=True)

BLUE = "#1a5fa8"
RED = "#c0392b"
GREEN = "#27ae60"
GRAY = "#555555"

nodes = [
    (0.0, 0.0),   # 0 corner
    (1.0, 0.0),   # 1 corner
    (1.0, 1.0),   # 2 corner
    (0.0, 1.0),   # 3 corner
    (0.5, 0.0),   # 4 edge midpoint
    (1.0, 0.5),   # 5 edge midpoint
    (0.5, 1.0),   # 6 edge midpoint
    (0.0, 0.5),   # 7 edge midpoint
    (0.5, 0.5),   # 8 center
]

# Colors: corners blue, edge midpoints red, center green
colors = [BLUE, BLUE, BLUE, BLUE, RED, RED, RED, RED, GREEN]

label_offsets = [
    (-0.09, -0.09),  # 0 bottom-left
    ( 0.06, -0.09),  # 1 bottom-right
    ( 0.06,  0.05),  # 2 top-right
    (-0.09,  0.05),  # 3 top-left
    ( 0.00, -0.10),  # 4 bottom edge
    ( 0.09,  0.00),  # 5 right edge
    ( 0.00,  0.08),  # 6 top edge
    (-0.09,  0.00),  # 7 left edge
    ( 0.09,  0.00),  # 8 center
]

coord_labels = [
    r"$(0,0)$",
    r"$(1,0)$",
    r"$(1,1)$",
    r"$(0,1)$",
    r"$(\frac{1}{2},0)$",
    r"$(1,\frac{1}{2})$",
    r"$(\frac{1}{2},1)$",
    r"$(0,\frac{1}{2})$",
    r"$(\frac{1}{2},\frac{1}{2})$",
]

fig, ax = plt.subplots(figsize=(4.8, 4.5))

# Square outline
sq = plt.Polygon([[0, 0], [1, 0], [1, 1], [0, 1]], closed=True,
                 edgecolor=GRAY, facecolor="#f0f4fa", linewidth=1.8)
ax.add_patch(sq)

for i, ((nx, ny), (ox, oy), clbl, color) in enumerate(
        zip(nodes, label_offsets, coord_labels, colors)):
    ax.plot(nx, ny, "o", color=color, markersize=8, zorder=5)
    ax.text(nx + ox, ny + oy + 0.09,
            f"$\\hat{{x}}_{i}$", fontsize=11, color=color,
            ha="center", va="bottom")
    ax.text(nx + ox, ny + oy - 0.01,
            clbl, fontsize=8, color=GRAY, ha="center", va="top")

ax.plot([], [], "o", color=BLUE, label="corner nodes")
ax.plot([], [], "o", color=RED, label="edge-midpoint nodes")
ax.plot([], [], "o", color=GREEN, label="center node")
ax.legend(fontsize=9, loc="upper right", frameon=False)

ax.set_xlim(-0.20, 1.20)
ax.set_ylim(-0.20, 1.22)
ax.set_aspect("equal")
ax.axis("off")

fig.tight_layout()
fig.savefig(out_path, format="svg", bbox_inches="tight")
print(f"Saved {out_path}")
