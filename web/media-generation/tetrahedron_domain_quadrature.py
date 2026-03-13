"""Generate tetrahedron-domain-quadrature.svg.

Shows quadrature points inside the reference tetrahedron for orders 1, 2, 3
(1, 4, 5 points). Each panel is a 3D wireframe view from the same viewpoint.
Positive weights are blue, negative are red. Marker size is proportional to
|weight|.
"""

from pathlib import Path
from itertools import combinations

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from yggdrasil.refdomains.tetrahedron import TetrahedronDomain

out_path = Path(__file__).parent.parent / "media" / "tetrahedron-domain-quadrature.svg"
out_path.parent.mkdir(parents=True, exist_ok=True)

BLUE = "#1a5fa8"
RED = "#c0392b"
ELEV, AZIM = 10, -100

# Reference tetrahedron vertices
VERTS = np.array([
    [0, 0, 0],
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1],
], dtype=float)

# Four triangular faces (vertex index triples)
FACES = [
    [VERTS[0], VERTS[1], VERTS[2]],
    [VERTS[0], VERTS[1], VERTS[3]],
    [VERTS[0], VERTS[2], VERTS[3]],
    [VERTS[1], VERTS[2], VERTS[3]],
]

domain = TetrahedronDomain()
orders = [1, 2, 3]
n_labels = {1: "1 pt", 2: "4 pts", 3: "5 pts"}

max_w = 1.0 / 6.0  # volume of reference tetrahedron

fig = plt.figure(figsize=(9, 3.4))

for idx, order in enumerate(orders, 1):
    ax = fig.add_subplot(1, 3, idx, projection="3d")
    ax.view_init(elev=ELEV, azim=AZIM)

    # Shaded faces (no edges here — drawn separately for clarity)
    faces = Poly3DCollection(FACES, alpha=0.07, facecolor="#1a5fa8",
                             edgecolor="none")
    ax.add_collection3d(faces)

    # Edges drawn explicitly for crispness
    for i, j in combinations(range(4), 2):
        xs, ys, zs = zip(VERTS[i], VERTS[j])
        ax.plot(xs, ys, zs, color="#333333", linewidth=1.6, zorder=2)

    # Vertices
    ax.scatter(*VERTS.T, color="#333333", s=22, zorder=5)

    # Quadrature points
    pts, wts = domain.quadrature(order)
    colors = [BLUE if w >= 0 else RED for w in wts]
    sizes = (np.abs(wts) / max_w) * 120
    ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2],
               s=sizes, c=colors, zorder=10, depthshade=False)

    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_zlim(0, 1)
    ax.set_xlabel("x", fontsize=8, labelpad=-4)
    ax.set_ylabel("y", fontsize=8, labelpad=-4)
    ax.set_zlabel("z", fontsize=8, labelpad=-4)
    ax.tick_params(labelsize=6, pad=-2)
    ax.set_title(f"order {order}  ({n_labels[order]})", fontsize=10)

fig.suptitle("Quadrature points on the reference tetrahedron", fontsize=11)
fig.tight_layout()
fig.savefig(out_path, format="svg", bbox_inches="tight")
print(f"Saved {out_path}")
