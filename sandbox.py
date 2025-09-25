import warnings
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

import PyOCN

G = nx.DiGraph()
n = 6
rows, cols = np.mgrid[0:n, 0:n]
rows, cols = rows[::-1].flatten(), cols.flatten()
for r, c in zip(rows, cols):
    G.add_node(r*n + c, pos=(r, c))
for r, c in zip(rows, cols):
    if r < n - 1:
        G.add_edge(r*n + c, (r + 1)*n + c)
    elif c < n - 1:
        G.add_edge(r*n + c, r*n + (c + 1))

ocn = PyOCN.OCN(init_structure=G, gamma=1, random_state=8472)

for _ in range(1000):
    ocn.single_erosion_event(temperature=0)
    PyOCN.validate_streamgraph(ocn.sg)
    # print(ocn.energy)
PyOCN.plot_streamgraph(ocn)
plt.show()



