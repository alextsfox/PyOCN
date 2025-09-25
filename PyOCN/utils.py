import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib.patches import ArrowStyle

from .ocn import OCN
from .streamgraph import StreamGraph

def plot_streamgraph(sg:StreamGraph|OCN, ax=None):
    """Plot the StreamGraph in matplotlib"""
    if isinstance(sg, OCN):
        sg = sg.sg
    return plot_positional_digraph(sg.dag, ax=ax)

def plot_positional_digraph(G:nx.DiGraph, ax=None):
    pos = nx.get_node_attributes(G, 'pos')
    # Transpose from (row, col) to (x, y) for plotting
    nrows = max(r for r, _ in pos.values()) + 1
    for node, (r, c) in pos.items():
        pos[node] = (c, nrows - r - 1)  

    drained_area = np.asarray(list(nx.get_node_attributes(G, 'drained_area').values()))
    drained_area = (drained_area - drained_area.min()) / (drained_area.max() - drained_area.min())
    size = 50 + drained_area*1000
    lw = 1 + drained_area*5
    if ax is None:
        _, ax = plt.subplots()
    nx.draw_networkx(G, pos=pos, node_size=size, width=lw, with_labels=False, arrowstyle=ArrowStyle("-|>", head_length=1.5, head_width=0.5), ax=ax)
    return ax

__all__ = [
    "plot_streamgraph",
    "plot_positional_digraph",
]