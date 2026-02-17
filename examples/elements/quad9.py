"""Plot the nine shape functions of the Quad9 element as 2D colormaps."""

import numpy as np
import matplotlib.pyplot as plt

from yggdrasil.elements import Quad9

element = Quad9()

n = 200
x1d = np.linspace(0.0, 1.0, n)
X, Y = np.meshgrid(x1d, x1d)
xi = np.column_stack([X.ravel(), Y.ravel()])
N = element.shape_functions(xi)
node_coords = element.node_coords

for i in range(element.num_nodes):
    fig, ax = plt.subplots(figsize=(5, 4.5))
    Z = N[:, i].reshape(n, n)
    vmin, vmax = Z.min(), Z.max()
    pc = ax.pcolormesh(X, Y, Z, cmap="viridis", vmin=vmin, vmax=vmax, shading="gouraud")
    levels = np.linspace(vmin, vmax, 11)
    ax.contour(X, Y, Z, levels=levels, colors="k", linewidths=0.5, alpha=0.4)
    fig.colorbar(pc, ax=ax)
    ax.plot(*node_coords.T, "ko", markersize=6)
    ax.set_xlabel("ξ")
    ax.set_ylabel("η")
    ax.set_title(f"Quad9 — Shape function N{i}")
    ax.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(f"quad9_n{i}.png", dpi=150)
    plt.close(fig)
    print(f"Wrote quad9_n{i}.png")
