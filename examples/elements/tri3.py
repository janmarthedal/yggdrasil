"""Plot the three shape functions of the Tri3 element as 2D colormaps."""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

from yggdrasil.elements import Tri3

element = Tri3()

# Create a fine triangulation over the reference triangle (0,0)-(1,0)-(0,1)
n = 50
xi_1d = np.linspace(0.0, 1.0, n)
x, y = [], []
for i in range(n):
    for j in range(n - i):
        x.append(xi_1d[i])
        y.append(xi_1d[j])
x = np.array(x)
y = np.array(y)
triang = mtri.Triangulation(x, y)

xi = np.column_stack([x, y])
N = element.shape_functions(xi)
node_coords = element.node_coords

for i in range(element.num_nodes):
    fig, ax = plt.subplots(figsize=(5, 4.5))
    tc = ax.tricontourf(triang, N[:, i], levels=np.linspace(0.0, 1.0, 21), cmap="viridis")
    ax.tricontour(triang, N[:, i], levels=np.linspace(0.0, 1.0, 11), colors="k",
                  linewidths=0.5, alpha=0.4)
    fig.colorbar(tc, ax=ax)
    ax.plot(*node_coords.T, "ko", markersize=6)
    ax.set_xlabel("ξ")
    ax.set_ylabel("η")
    ax.set_title(f"Tri3 — Shape function N{i}")
    ax.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(f"tri3_n{i}.png", dpi=150)
    plt.close(fig)
    print(f"Wrote tri3_n{i}.png")
