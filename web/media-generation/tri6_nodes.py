"""Generate tri6-nodes.svg.

Shows the reference triangle with all 6 nodes of the Tri6 element
labelled by index and coordinate.
"""

from pathlib import Path

import matplotlib.pyplot as plt

out_path = Path(__file__).parent.parent / "media" / "tri6-nodes.svg"
out_path.parent.mkdir(parents=True, exist_ok=True)

BLUE = "#1a5fa8"
GRAY = "#555555"

nodes = [
    (0.0, 0.0),
    (1.0, 0.0),
    (0.0, 1.0),
    (0.5, 0.0),
    (0.5, 0.5),
    (0.0, 0.5),
]

# (dx, dy) nudge for each label so they don't overlap the marker
label_offsets = [
    (-0.07, -0.08),  # 0: bottom-left corner
    ( 0.05, -0.08),  # 1: bottom-right corner
    (-0.07,  0.05),  # 2: top corner
    ( 0.00, -0.09),  # 3: bottom edge midpoint
    ( 0.06,  0.02),  # 4: hypotenuse midpoint
    (-0.07,  0.02),  # 5: left edge midpoint
]

coord_labels = [
    r"$(0,0)$",
    r"$(1,0)$",
    r"$(0,1)$",
    r"$(\frac{1}{2},0)$",
    r"$(\frac{1}{2},\frac{1}{2})$",
    r"$(0,\frac{1}{2})$",
]

fig, ax = plt.subplots(figsize=(4.5, 4.2))

# Triangle outline
tri = plt.Polygon([[0, 0], [1, 0], [0, 1]], closed=True,
                  edgecolor=GRAY, facecolor="#f0f4fa", linewidth=1.8)
ax.add_patch(tri)

for i, ((nx, ny), (ox, oy), clbl) in enumerate(
        zip(nodes, label_offsets, coord_labels)):
    color = BLUE if i < 3 else "#c0392b"
    ax.plot(nx, ny, "o", color=color, markersize=8, zorder=5)
    ax.text(nx + ox, ny + oy + 0.09,
            f"$\\hat{{x}}_{i}$", fontsize=11, color=color,
            ha="center", va="bottom")
    ax.text(nx + ox, ny + oy - 0.01,
            clbl, fontsize=8, color=GRAY, ha="center", va="top")

# Legend
ax.plot([], [], "o", color=BLUE, label="vertex nodes")
ax.plot([], [], "o", color="#c0392b", label="edge-midpoint nodes")
ax.legend(fontsize=9, loc="upper right", frameon=False)

ax.set_xlim(-0.18, 1.12)
ax.set_ylim(-0.18, 1.18)
ax.set_aspect("equal")
ax.axis("off")

fig.tight_layout()
fig.savefig(out_path, format="svg", bbox_inches="tight")
print(f"Saved {out_path}")
