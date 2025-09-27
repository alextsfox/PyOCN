from __future__ import annotations

from typing import Any, TYPE_CHECKING

import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import warnings

if TYPE_CHECKING:
    from .ocn import OCN

def _pos_to_xy(dag:nx.DiGraph) -> dict[Any, tuple[float, float]]:
    """Convert 'pos' attributes from (row, col) to (x, y) for plotting."""
    pos = nx.get_node_attributes(dag, 'pos')
    nrows = max(r for r, _ in pos.values()) + 1
    for node, (r, c) in pos.items():
        pos[node] = (c, nrows - r - 1)  
    return pos

def plot_ocn_as_dag(ocn:OCN, attribute=None, ax=None, norm=None, **kwargs):
    """Plot the OCN as a directed acyclic graph (DAG) using networkx.
    
    Parameters
    ----------
    ocn: OCN
        The OCN instance to plot.
    attribute: str, optional
        Node attribute to coloring nodes by (e.g., 'drained_area' or 'energy').
    ax: matplotlib axes, optional
        Axes to plot on. If None, a new figure and axes are created.
    norm: matplotlib.colors.Normalize, optional
        Normalization for node colors if attribute is specified.
    **kwargs: additional keyword arguments (e.g. cmap, vmin, vmax).
    """
    
    dag = ocn.to_digraph()
    pos = _pos_to_xy(dag)

    if ax is None:
        _, ax = plt.subplots()

    node_color = "C0"
    if attribute is not None:
        node_color = list(nx.get_node_attributes(dag, attribute).values())

    if norm is not None:
        if ("vmin" in kwargs or "vmax" in kwargs):
            warnings.warn("norm is specified, ignoring vmin/vmax.")
        kwargs["vmin"] = 0
        kwargs["vmax"] = 1
        node_color = norm(node_color)

    p = nx.draw_networkx(dag, node_color=node_color, pos=pos, ax=ax, **kwargs)
    return p, ax

def plot_ocn_energy_raster(ocn:OCN, ax=None, **kwargs):
    """Plot a raster of the OCN, colored by cell energy.
    
    Parameters
    ----------
    ocn: OCN
        The OCN instance to plot.
    ax: matplotlib axes, optional
        Axes to plot on. If None, a new figure and axes are created.
    **kwargs: additional keyword arguments passed to pcolormesh.
    """

    dag = ocn.to_digraph()
    energy = np.zeros(ocn.dims)
    for node in dag.nodes:
        r, c = dag.nodes[node]['pos']
        energy[r, c] = dag.nodes[node]['energy']

    if "cmap" not in kwargs:
        kwargs["cmap"] = "terrain"
        
    if ax is None:
        _, ax = plt.subplots()

    ax.imshow(energy, **kwargs)
    return ax
    

def plot_positional_digraph(dag:nx.DiGraph, ax=None, **kwargs):
    """Plot a directed acyclic graph (DAG) with nodes positioned according to their 'pos' attribute and sized by drained area."""
    pos = _pos_to_xy(dag)

    if ax is None:
        _, ax = plt.subplots()

    p = nx.draw_networkx(dag, pos=pos, ax=ax, **kwargs)
    return p, ax