"""
Plotting utilities for PyOCN.

This module provides convenience functions to visualize OCNs using Matplotlib
and NetworkX.

Notes
-----
- Node positions are stored as ``pos=(row, col)`` in the data model. For
    visualization, these are converted to cartesian coordinates so that 
    (row, col) = (0, 0) is at the bottom-left corner.
"""

from __future__ import annotations

from typing import Any, TYPE_CHECKING

import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import warnings

if TYPE_CHECKING:
    from .ocn import OCN

def _pos_to_xy(dag: nx.DiGraph) -> dict[Any, tuple[float, float]]:
    """
    Convert node ``pos`` from (row, col) to plotting coordinates (x, y).

    Parameters
    ----------
    dag : nx.DiGraph
        Graph whose nodes have ``pos=(row, col)`` attributes.

    Returns
    -------
    dict[Any, tuple[float, float]]
        Mapping from node to ``(x, y)`` where ``x=col`` and ``y`` is flipped
        to plot with origin at the bottom.

    Notes
    -----
    The vertical coordinate is transformed as ``y = nrows - row - 1`` to
    match typical plotting conventions.
    """
    pos = nx.get_node_attributes(dag, 'pos')
    nrows = max(r for r, _ in pos.values()) + 1
    for node, (r, c) in pos.items():
        pos[node] = (c, nrows - r - 1)  
    return pos

def plot_ocn_as_dag(ocn: OCN, attribute: str | None = None, ax=None, norm=None, **kwargs):
    """
    Plot the OCN as a DAG using NetworkX.

    Parameters
    ----------
    ocn : OCN
        The OCN instance to plot.
    attribute : str, optional
        Node attribute name for coloring (e.g., ``'drained_area'`` or
        ``'energy'``). If omitted, a uniform color is used.
    ax : matplotlib.axes.Axes, optional
        Target axes. If ``None``, a new figure and axes are created.
    norm : matplotlib.colors.Normalize, optional
        Normalization to apply to node colors when ``attribute`` is provided.
        If specified, any ``vmin``/``vmax`` passed in ``**kwargs`` are ignored.
    **kwargs
        Additional keyword arguments forwarded to
        :func:`networkx.draw_networkx` (e.g., ``cmap``, ``vmin``, ``vmax``,
        size and style options).

    Returns
    -------
    tuple
        A pair ``(artists, ax)`` where ``artists`` is the object returned by
        ``networkx.draw_networkx`` (backend-dependent; often ``None``) and
        ``ax`` is the axes used for drawing.
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

def plot_ocn_raster(ocn: OCN, attribute:str='energy', ax=None, **kwargs):
    """
    Plot a raster image of grid cell energies.

    Parameters
    ----------
    ocn : OCN
        The OCN instance to plot.
    attribute : str, default 'energy'
        The node attribute to visualize (e.g., 'energy', 'drained_area').
    ax : matplotlib.axes.Axes, optional
        Target axes. If ``None``, a new figure and axes are created.
    **kwargs
        Additional keyword arguments forwarded to ``imshow`` (e.g.,
        ``cmap``, interpolation options, and normalization).

    Returns
    -------
    matplotlib.axes.Axes
        The axes containing the rendered image.
    """

    dag = ocn.to_digraph()
    energy = np.zeros(ocn.dims)
    for node in dag.nodes:
        r, c = dag.nodes[node]['pos']
        energy[r, c] = dag.nodes[node][attribute]

    if "cmap" not in kwargs:
        kwargs["cmap"] = "terrain_r"
        
    if ax is None:
        _, ax = plt.subplots()

    ax.imshow(energy, **kwargs)
    return ax
    

def plot_positional_digraph(dag: nx.DiGraph, ax=None, **kwargs):
    """
    Plot a DAG with node positions taken from their ``pos`` attributes.

    Parameters
    ----------
    dag : nx.DiGraph
        Graph whose nodes have ``pos=(row, col)``.
    ax : matplotlib.axes.Axes, optional
        Target axes. If ``None``, a new figure and axes are created.
    **kwargs
        Additional keyword arguments forwarded to
        :func:`networkx.draw_networkx`.

    Returns
    -------
    tuple
        A pair ``(artists, ax)`` where ``artists`` is the object returned by
        ``networkx.draw_networkx`` (often ``None``) and ``ax`` is the axes
        used for drawing.
    """
    pos = _pos_to_xy(dag)

    if ax is None:
        _, ax = plt.subplots()

    p = nx.draw_networkx(dag, pos=pos, ax=ax, **kwargs)
    return p, ax