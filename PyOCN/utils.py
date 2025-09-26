import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib.patches import ArrowStyle

from .ocn import OCN

def plot_ocn_as_dag(ocn:OCN, ax=None):
    dag = ocn.to_digraph()
    return plot_positional_digraph(dag, ax=ax)

def plot_positional_digraph(dag:nx.DiGraph, ax=None):
    pos = nx.get_node_attributes(dag, 'pos')

    # Transpose from (row, col) to (x, y) for plotting
    nrows = max(r for r, _ in pos.values()) + 1
    for node, (r, c) in pos.items():
        pos[node] = (c, nrows - r - 1)  

    drained_areas = dict()
    for node in nx.topological_sort(dag):
        succs = list(dag.successors(node))
        drained_areas[node] = 1 + sum(drained_areas[s] for s in succs)
    
    sizes = {k: 50 + da*1000 for k, da in drained_areas.items()}
    lws = {k: 1 + da*5 for k, da in drained_areas.items()}
    if ax is None:
        _, ax = plt.subplots()
    nx.draw_networkx(dag, pos=pos, node_size=list(sizes.values()), width=list(lws.values()), with_labels=False, arrowstyle=ArrowStyle("-|>", head_length=1.5, head_width=0.5), ax=ax)
    return ax