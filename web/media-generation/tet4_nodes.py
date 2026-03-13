"""Generate tet4-nodes.svg.

Shows the reference tetrahedron with all 4 nodes of the Tet4 element
labelled by index and coordinate.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

out_path = Path(__file__).parent.parent / "media" / "tet4-nodes.svg"
out_path.parent.mkdir(parents=True, exist_ok=True)

BLUE = "#1a5fa8"
GRAY = "#555555"
FACE_COLOR = "#d0dff0"
EDGE_COLOR = "#555555"

nodes = np.array([
    [0.0, 0.0, 0.0],  # 0
    [1.0, 0.0, 0.0],  # 1
    [0.0, 1.0, 0.0],  # 2
    [0.0, 0.0, 1.0],  # 3
])

node_labels = [
    (r"$\hat{x}_0$", r"$(0,0,0)$"),
    (r"$\hat{x}_1$", r"$(1,0,0)$"),
    (r"$\hat{x}_2$", r"$(0,1,0)$"),
    (r"$\hat{x}_3$", r"$(0,0,1)$"),
]

# Nudge directions to avoid overlap with the mesh
label_offsets = [
    (-0.10, -0.10, -0.10),  # 0: origin — push toward negative axes
    ( 0.08,  0.00, -0.08),  # 1: x-axis tip
    (-0.10,  0.08, -0.08),  # 2: y-axis tip
    (-0.10, -0.08,  0.08),  # 3: z-axis tip
]

# Four triangular faces (each defined by 3 node indices)
faces = [
    [0, 1, 2],  # base (z=0)
    [0, 1, 3],  # front (y=0)
    [0, 2, 3],  # left (x=0)
    [1, 2, 3],  # hypotenuse face
]

fig = plt.figure(figsize=(4.5, 4.5))
ax = fig.add_subplot(111, projection="3d")

# Draw faces
polys = [[nodes[i] for i in face] for face in faces]
face_alphas = [0.25, 0.20, 0.20, 0.30]
for poly, alpha in zip(polys, face_alphas):
    col = Poly3DCollection([poly], alpha=alpha,
                           facecolor=FACE_COLOR, edgecolor=EDGE_COLOR,
                           linewidth=1.2)
    ax.add_collection3d(col)

# Draw nodes
ax.scatter(*nodes.T, color=BLUE, s=60, zorder=5)

# Labels
for i, ((nx, ny, nz), (ox, oy, oz), (name, coord)) in enumerate(
        zip(nodes, label_offsets, node_labels)):
    ax.text(nx + ox, ny + oy, nz + oz + 0.08,
            name, fontsize=11, color=BLUE, ha="center", va="bottom")
    ax.text(nx + ox, ny + oy, nz + oz - 0.01,
            coord, fontsize=8, color=GRAY, ha="center", va="top")

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
