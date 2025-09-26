import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.patches import ArrowStyle

from .ocn import OCN

def plot_ocn_as_dag(ocn:OCN, ax=None):
    """Plot the OCN as a directed acyclic graph (DAG) using networkx."""
    dag = ocn.to_digraph()
    return plot_positional_digraph(dag, ax=ax)

def plot_positional_digraph(dag:nx.DiGraph, ax=None):
    """Plot a directed acyclic graph (DAG) with nodes positioned according to their 'pos' attribute and sized by drained area."""
    pos = nx.get_node_attributes(dag, 'pos')

    # Transpose from (row, col) to (x, y) for plotting
    nrows = max(r for r, _ in pos.values()) + 1
    for node, (r, c) in pos.items():
        pos[node] = (c, nrows - r - 1)  

    drained_areas = dict()
    for node in nx.topological_sort(dag):
        preds = list(dag.predecessors(node))
        drained_areas[node] = 1 + sum(drained_areas[p] for p in preds)
    max_area = max(drained_areas.values())
    drained_areas = {k: v/max_area for k, v in drained_areas.items()}

    sizes = {k: 50 + da*1000 for k, da in drained_areas.items()}
    lws = {k: 1 + da*5 for k, da in drained_areas.items()}
    if ax is None:
        _, ax = plt.subplots()
    
    p = nx.draw_networkx(
        dag, 
        pos=pos, 
        node_size=list(sizes.values()), 
        width=list(lws.values()), 
        with_labels=False, 
        arrowstyle=ArrowStyle("-|>", head_length=1.5, head_width=0.5), 
        ax=ax
    )
    return p, ax