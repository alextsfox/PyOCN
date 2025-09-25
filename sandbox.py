import warnings
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

import PyOCN

################################################
# basic 4x4 graph, no issues
label = "basic 4x4"
G = nx.DiGraph()

G.add_node(0, pos=(0, 0), drained_area=1)
G.add_node(1, pos=(0, 1), drained_area=1)
G.add_node(2, pos=(0, 2), drained_area=1)
G.add_node(3, pos=(0, 3), drained_area=1)
G.add_node(4, pos=(1, 0), drained_area=2)
G.add_node(5, pos=(1, 1), drained_area=2)
G.add_node(6, pos=(1, 2), drained_area=2)
G.add_node(7, pos=(1, 3), drained_area=2)
G.add_node(8, pos=(2, 0), drained_area=3)
G.add_node(9, pos=(2, 1), drained_area=3)
G.add_node(10, pos=(2, 2), drained_area=3)
G.add_node(11, pos=(2, 3), drained_area=3)
G.add_node(12, pos=(3, 0), drained_area=4)
G.add_node(13, pos=(3, 1), drained_area=8)
G.add_node(14, pos=(3, 2), drained_area=12)
G.add_node(15, pos=(3, 3), drained_area=16)

G.add_edge(0, 4)
G.add_edge(1, 5)
G.add_edge(2, 6)
G.add_edge(3, 7)
G.add_edge(4, 8)
G.add_edge(5, 9)
G.add_edge(6, 10)
G.add_edge(7, 11)
G.add_edge(8, 12)
G.add_edge(9, 13)
G.add_edge(10, 14)
G.add_edge(11, 15)
G.add_edge(12, 13)
G.add_edge(13, 14)
G.add_edge(14, 15)

ocn = PyOCN.OCN(init_structure=G, gamma=1, random_state=8472)

PyOCN.plot_streamgraph(ocn)
plt.show()