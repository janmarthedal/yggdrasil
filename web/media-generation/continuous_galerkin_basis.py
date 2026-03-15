"""Generate continuous-galerkin-basis.svg.

Shows a piecewise-linear basis function φ_i on a triangular mesh whose outer
boundary nodes lie on a square but whose interior nodes are placed irregularly.
The elements sharing the featured interior node i are shaded by function value
(1 at node i, 0 at all neighbouring nodes). The patch boundary is highlighted
in blue to illustrate continuity enforcement through shared nodes.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.collections as mc
import matplotlib.tri as mtri
import numpy as np
from scipy.spatial import Delaunay

out_path = Path(__file__).parent.parent / "media" / "continuous-galerkin-basis.svg"
out_path.parent.mkdir(parents=True, exist_ok=True)

BLUE = "#1a5fa8"
GRAY = "#555555"
OUTER_FILL = "#e4e4e4"

# ---------------------------------------------------------------------------
# Nodes: boundary nodes lie exactly on [0,4]², interior nodes are irregular.
# The featured node is the first interior node (near the centre).
# ---------------------------------------------------------------------------
boundary_nodes = np.array([
    # corners
    [0.0, 0.0], [4.0, 0.0], [4.0, 4.0], [0.0, 4.0],
    # bottom edge (y=0)
    [1.2, 0.0], [2.8, 0.0],
    # right edge (x=4)
    [4.0, 1.5], [4.0, 3.0],
    # top edge (y=4)
    [2.6, 4.0], [1.0, 4.0],
    # left edge (x=0)
    [0.0, 2.7], [0.0, 1.1],
])

interior_nodes = np.array([
    [2.1, 2.0],  # 12 — featured node (near centre)
    [1.0, 1.1],  # 13
    [3.0, 1.2],  # 14
    [3.3, 2.9],  # 15
    [1.1, 3.2],  # 16
    [2.2, 3.3],  # 17
])

nodes = np.vstack([boundary_nodes, interior_nodes])
featured = len(boundary_nodes)  # index of the featured node

# ---------------------------------------------------------------------------
# Delaunay triangulation
# ---------------------------------------------------------------------------
tri = Delaunay(nodes)
all_tris = tri.simplices

# Split into support (contain featured node) and outer
support_mask = np.any(all_tris == featured, axis=1)
support_tris = all_tris[support_mask]
outer_tris = all_tris[~support_mask]

# ---------------------------------------------------------------------------
# Patch-boundary edges (shared between support and outer triangles)
# ---------------------------------------------------------------------------
def edge_set(tris):
    edges = set()
    for t in tris:
        for k in range(3):
            edges.add(tuple(sorted([t[k], t[(k + 1) % 3]])))
    return edges

patch_boundary = edge_set(support_tris) & edge_set(outer_tris)

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
x, y = nodes[:, 0], nodes[:, 1]

fig, ax = plt.subplots(figsize=(5.0, 5.0))

# 1. Outer triangles — light gray fill
outer_col = mc.PolyCollection(nodes[outer_tris], facecolor=OUTER_FILL,
                               edgecolor="none", zorder=1)
ax.add_collection(outer_col)

# 2. Support triangles — coloured by φ value.
# tricontourf with many levels approximates a smooth gradient in SVG
# (tripcolor gouraud shading is not supported by the SVG backend).
phi = np.zeros(len(nodes))
phi[featured] = 1.0
support_triang = mtri.Triangulation(x, y, support_tris)
levels = np.linspace(0.0, 1.0, 50)
tcf = ax.tricontourf(support_triang, phi, levels=levels, cmap="Blues",
                     vmin=0.0, vmax=1.0, zorder=2)

# 3. All mesh edges — thin gray
drawn = set()
for t in all_tris:
    for k in range(3):
        edge = tuple(sorted([t[k], t[(k + 1) % 3]]))
        if edge not in drawn:
            a, b = edge
            ax.plot([x[a], x[b]], [y[a], y[b]], color=GRAY,
                    linewidth=0.9, zorder=3)
            drawn.add(edge)

# 4. Patch boundary — highlighted in blue
for a, b in patch_boundary:
    ax.plot([x[a], x[b]], [y[a], y[b]], color=BLUE, linewidth=2.2, zorder=4)

# 5. Non-featured nodes
others = [i for i in range(len(nodes)) if i != featured]
ax.plot(x[others], y[others], "o", color=GRAY, markersize=4, zorder=5)

# 6. Featured node
ax.plot(x[featured], y[featured], "o", color=BLUE, markersize=9, zorder=6)
ax.text(x[featured] + 0.1, y[featured] + 0.15, r"$x_i$", fontsize=13,
        color=BLUE, zorder=7)

# 7. Colorbar
cbar = fig.colorbar(tcf, ax=ax, fraction=0.04, pad=0.02, aspect=20)
cbar.set_label(r"$\varphi_i$", fontsize=13)
cbar.set_ticks([0, 0.5, 1])

ax.set_aspect("equal")
ax.axis("off")
fig.tight_layout()
fig.savefig(out_path, format="svg", bbox_inches="tight")
print(f"Saved {out_path}")
