from itertools import product

from tqdm import trange, tqdm
import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx

import PyOCN

################################################
# basic 4x4 graph, no issues
# label = "basic 4x4"
# G = nx.DiGraph()

# rows, cols = 8, 12

# # H shape shape
# halfway = cols // 2
# for i, j in product(range(rows), range(cols)):
#     n = i*cols + j
#     G.add_node(n, pos=(i, j))
#     if i == j and i > 0:  # main diagonal
#         G.add_edge(n, n - cols - 1)
#     elif i > j:
#         G.add_edge(n, n - cols)
#     elif j > i:
#         G.add_edge(n, n - 1)


ocn = PyOCN.OCN.from_net_type(
    "H", 
    dims=(26, 26), 
    gamma=0.1, 
    annealing_schedule=1.0, 
    random_state=8473
)

for _ in trange(30000):
    ocn.single_erosion_event(0.999)

dag = ocn.to_digraph()
print(PyOCN._streamgraph_convert.validate_digraph(dag))

PyOCN.plot_ocn_energy_raster(ocn=ocn)
# PyOCN.plot_ocn_as_dag(ocn, attribute='energy', with_labels=False, node_size=10)
plt.show()