"""Plot the three shape functions of the Line3 element, one SVG per node."""

import numpy as np
import matplotlib.pyplot as plt

from yggdrasil.elements import Line3

element = Line3()
xi = np.linspace(0.0, 1.0, 200).reshape(-1, 1)
N = element.shape_functions(xi)
node_coords = element.node_coords[:, 0]

for i in range(element.num_nodes):
    fig, ax = plt.subplots(figsize=(5, 3))
    ax.plot(xi[:, 0], N[:, i], linewidth=2)
    ax.scatter(node_coords, [1.0 if j == i else 0.0 for j in range(element.num_nodes)],
               color="black", zorder=5)
    ax.set_xlabel("ξ")
    ax.set_ylabel(f"N{i}(ξ)")
    ax.set_title(f"Line3 — Shape function N{i}")
    ax.set_ylim(-0.15, 1.15)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"line3_n{i}.svg")
    plt.close(fig)
    print(f"Wrote line3_n{i}.svg")
