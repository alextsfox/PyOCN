from itertools import product

import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx

import PyOCN

################################################
# basic 4x4 graph, no issues
label = "basic 4x4"
G = nx.DiGraph()

rows, cols = 8, 12

# H shape shape
halfway = cols // 2
for i, j in product(range(rows), range(cols)):
    n = i*cols + j
    G.add_node(n, pos=(i, j))
    if i == j and i > 0:  # main diagonal
        G.add_edge(n, n - cols - 1)
    elif i > j:
        G.add_edge(n, n - cols)
    elif j > i:
        G.add_edge(n, n - 1)


ocn = PyOCN.OCN.from_digraph(G, random_state=8472)
PyOCN.plot_ocn_as_dag(ocn, attribute='energy', cmap='Blues')
plt.show()

