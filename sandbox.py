from itertools import product

import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx

import PyOCN

################################################
# basic 4x4 graph, no issues
label = "basic 4x4"
G = nx.DiGraph()

rows, cols = 10, 12

# I shape
halfway = cols // 2
for i, j in product(range(rows), range(cols)):
    n = i*cols + j
    G.add_node(n, pos=(i, j))

    if j < halfway:
        G.add_edge(n, n+1)
    elif j > halfway:
        G.add_edge(n, n-1)
    elif i < rows - 1:
        G.add_edge(n, n + cols)

ocn = PyOCN.OCN.from_digraph(G, random_state=8472)
PyOCN.plot_ocn_as_dag(ocn, attribute='energy', cmap='Blues')
plt.show()