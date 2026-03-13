"""Generate hexahedron-domain-quadrature.svg.

Shows the tensor-product Gauss-Legendre quadrature points on [0,1]^3
for order 5 (3x3x3 = 27 points). Marker size is proportional to weight.
"""

from pathlib import Path
from itertools import combinations

import matplotlib.pyplot as plt
import numpy as np

from yggdrasil.refdomains.hexahedron import HexahedronDomain

out_path = Path(__file__).parent.parent / "media" / "hexahedron-domain-quadrature.svg"
out_path.parent.mkdir(parents=True, exist_ok=True)

BLUE = "#1a5fa8"
ELEV, AZIM = 18, -105

# Reference cube vertices and edges
VERTS = np.array([[x, y, z] for x in [0, 1] for y in [0, 1] for z in [0, 1]])
EDGES = [(i, j) for i, j in combinations(range(8), 2)
         if sum(abs(VERTS[i] - VERTS[j])) == 1]

domain = HexahedronDomain()
pts, wts = domain.quadrature(5)
n = int(round(len(pts) ** (1/3)))

fig = plt.figure(figsize=(4.5, 4.5))
ax = fig.add_subplot(111, projection="3d")
ax.view_init(elev=ELEV, azim=AZIM)

# Draw edges
for i, j in EDGES:
    xs, ys, zs = zip(VERTS[i], VERTS[j])
    ax.plot(xs, ys, zs, color="#333333", linewidth=1.4, zorder=2)

# Draw vertices
ax.scatter(*VERTS.T, color="#333333", s=18, zorder=5)

# Draw quadrature points
sizes = (wts / wts.max()) * 120
ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2],
           s=sizes, color=BLUE, zorder=10, depthshade=False)

ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_zlim(0, 1)
ax.set_xlabel("x", fontsize=9, labelpad=-4)
ax.set_ylabel("y", fontsize=9, labelpad=-4)
ax.set_zlabel("z", fontsize=9, labelpad=-4)
ax.tick_params(labelsize=7, pad=-2)
ax.set_title(f"Order 5  ({n}×{n}×{n} = {n**3} points)", fontsize=11)

fig.tight_layout()
fig.savefig(out_path, format="svg", bbox_inches="tight")
print(f"Saved {out_path}")
