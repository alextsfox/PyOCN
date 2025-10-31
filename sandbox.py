import PyOCN as po
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

ocn = po.OCN.from_net_type("E", dims=(32, 32), random_state=8472)
ocn.fit(pbar=True)
G = ocn.to_digraph()
for n in G.nodes:
    if G.in_degree(n) > 2:
        print(n)
        break
n = 5
po.utils.get_subwatersheds(G, node=n)
# po.plotting.plot_ocn_as_dag(ocn, ax=ax, node_size=5, with_labels="False")