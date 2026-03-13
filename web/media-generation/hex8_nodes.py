"""Generate hex8-nodes.svg.

Shows the reference cube with all 8 nodes of the Hex8 element
labelled by index and coordinate.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

out_path = Path(__file__).parent.parent / "media" / "hex8-nodes.svg"
out_path.parent.mkdir(parents=True, exist_ok=True)

BLUE = "#1a5fa8"
RED = "#c0392b"
GRAY = "#555555"
FACE_COLOR = "#d0dff0"
EDGE_COLOR = "#555555"

nodes = np.array([
    [0.0, 0.0, 0.0],  # 0 bottom
    [1.0, 0.0, 0.0],  # 1 bottom
    [1.0, 1.0, 0.0],  # 2 bottom
    [0.0, 1.0, 0.0],  # 3 bottom
    [0.0, 0.0, 1.0],  # 4 top
    [1.0, 0.0, 1.0],  # 5 top
    [1.0, 1.0, 1.0],  # 6 top
    [0.0, 1.0, 1.0],  # 7 top
])

# Colors: bottom face blue, top face red
colors = [BLUE] * 4 + [RED] * 4

coord_labels = [
    r"$(0,0,0)$",
    r"$(1,0,0)$",
    r"$(1,1,0)$",
    r"$(0,1,0)$",
    r"$(0,0,1)$",
    r"$(1,0,1)$",
    r"$(1,1,1)$",
    r"$(0,1,1)$",
]

label_offsets = [
    (-0.10, -0.10, -0.10),  # 0
    ( 0.08, -0.10, -0.10),  # 1
    ( 0.08,  0.08, -0.10),  # 2
    (-0.10,  0.08, -0.10),  # 3
    (-0.10, -0.10,  0.08),  # 4
    ( 0.08, -0.10,  0.08),  # 5
    ( 0.08,  0.08,  0.08),  # 6
    (-0.10,  0.08,  0.08),  # 7
]

# Six faces of the cube
faces = [
    [0, 1, 2, 3],  # bottom
    [4, 5, 6, 7],  # top
    [0, 1, 5, 4],  # front
    [2, 3, 7, 6],  # back
    [0, 3, 7, 4],  # left
    [1, 2, 6, 5],  # right
]

fig = plt.figure(figsize=(5.0, 4.8))
ax = fig.add_subplot(111, projection="3d")

polys = [[nodes[i] for i in face] for face in faces]
col = Poly3DCollection(polys, alpha=0.12,
                       facecolor=FACE_COLOR, edgecolor=EDGE_COLOR,
                       linewidth=1.0)
ax.add_collection3d(col)

ax.scatter(*nodes.T, c=colors, s=60, zorder=5)

for i, ((nx, ny, nz), (ox, oy, oz), clbl, color) in enumerate(
        zip(nodes, label_offsets, coord_labels, colors)):
    ax.text(nx + ox, ny + oy, nz + oz + 0.08,
            f"$\\hat{{x}}_{i}$", fontsize=10, color=color,
            ha="center", va="bottom")
    ax.text(nx + ox, ny + oy, nz + oz - 0.01,
            clbl, fontsize=7, color=GRAY, ha="center", va="top")

ax.plot([], [], "o", color=BLUE, label="bottom face")
ax.plot([], [], "o", color=RED, label="top face")
ax.legend(fontsize=9, loc="upper left", frameon=False)

ax.set_xlim(-0.05, 1.05)
ax.set_ylim(-0.05, 1.05)
ax.set_zlim(-0.05, 1.05)
ax.set_xlabel(r"$\hat{x}$", labelpad=2)
ax.set_ylabel(r"$\hat{y}$", labelpad=2)
ax.set_zlabel(r"$\hat{z}$", labelpad=2)
ax.tick_params(labelsize=7)
ax.set_box_aspect([1, 1, 1])
ax.view_init(elev=22, azim=-55)

fig.tight_layout()
fig.savefig(out_path, format="svg", bbox_inches="tight")
print(f"Saved {out_path}")
