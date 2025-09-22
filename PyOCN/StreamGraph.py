"""
StreamGraph.py

High-level Python interface for the StreamGraph C library.

Author: Alexander S Fox
Copyright: (c) 2025 Alexander S Fox. All rights reserved.

This file is part of the OCN project.
"""

from ctypes import byref
from dataclasses import dataclass
import itertools
import ctypes

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.patches import ArrowStyle

from . import _StreamGraphC as _sgc

# """Python representation of a vertex in the StreamGraph....might not use"""
# @dataclass(slots=True)
# class Vertex:
#     drained_area: np.uint32
#     adown: np.uint32
#     edges: np.uint8
#     downstream: np.uint8
#     visited: np.uint8

#     @property
#     def is_root(self) -> bool:
#         return bool(self.downstream == _sgc.libocn.IS_ROOT)

class StreamGraph:
    def __init__(self, shape:tuple[int, int]=None, root:tuple[int, int]=None, init_structure="test"):
        """
        Initialize a StreamGraph object without any structure.

        Args:
            shape (tuple[int, int]): The dimensions of the graph as (rows, columns).
            root (tuple[int, int], optional): The (row, column) index of the root vertex. 
                If not provided, defaults to the bottom-right corner of the grid.

        Creates a directed graph of vertices arranged in a 2D grid, where each vertex has a drained area
        and flows downstream to another vertex. Initializes the underlying C StreamGraph structure and 
        prepares the vertex array for manipulation and analysis.
        """
        if init_structure == "test":
            self._c_graph = _sgc.libocn.sg_make_test_graph()

        else:    
            self._c_graph = _sgc.StreamGraphC()
            if shape is None:
                raise ValueError("Must provide shape when init_structure is not 'test'")
            i_root, j_root = shape[0]-1, shape[1]-1
            if root is not None:
                i_root, j_root = root

            status = _sgc.libocn.sg_create(byref(self._c_graph), shape[0], shape[1], i_root, j_root)
            if _sgc.STATUS_CODES[status] != "SUCCESS":
                raise RuntimeError(f"Failed to create StreamGraph: {_sgc.STATUS_CODES.get(status, 'Unknown error code')} ({status})")
        

    def __getitem__(self, idx):
        return self.vertices[idx]
    def __repr__(self):
        return repr(self._c_graph)
    def __str__(self):
        return f"StreamGraph(m={self.m}, n={self.n}, i_root={self.i_root}, j_root={self.j_root}, energy={self.energy:.5f}, vertices=<{self.m * self.n} vertices>)"
    def __del__(self):
        _sgc.libocn.sg_destroy(byref(self._c_graph))
    
    @property
    def m(self) -> np.uint16:
        return np.uint16(self._c_graph.m)
    @property
    def n(self) -> np.uint16:
        return np.uint16(self._c_graph.n)
    @property
    def shape(self) -> tuple[int, int]:
        return (self.m, self.n)
    @property
    def i_root(self) -> np.uint16:
        return np.uint16(self._c_graph.i_root)
    @property
    def j_root(self) -> np.uint16:
        return np.uint16(self._c_graph.j_root)
    @property
    def energy(self) -> np.float64:
        return np.float64(self._c_graph.energy)
    
    # TODO: figure out how to make this compatible with our in-memory tiling representation (not implemented currently in streamgraph.c)
    # TODO: vectorize without warnings?????
    @property
    def vertices(self) -> np.ndarray:
        # vertices = np.ctypeslib.as_array(self._c_graph.vertices, shape=(self.m*self.n,))
        vert_c = _sgc.VertexC()
        vertices = np.empty((self.m*self.n,), dtype=_sgc.vert_dtype)
        for a in range(self.m*self.n):
            _ = _sgc.libocn.sg_get_lin(self._c_graph, byref(vert_c), a)
            vertices[a] = (vert_c.drained_area, vert_c.adown, vert_c.edges, vert_c.downstream, vert_c.visited)
        return vertices
    
    def single_erosion_event(self):
        """Perform a single erosion event on the StreamGraph."""
        status = _sgc.libocn.sg_single_erosion_event(byref(self._c_graph))
        if _sgc.STATUS_CODES[status] != "SUCCESS":
            raise RuntimeError(f"Failed to perform single erosion event: {_sgc.STATUS_CODES.get(status, 'Unknown error code')} ({status})")
        
    def to_networkx(self) -> nx.DiGraph:
        """Convert the StreamGraph to a NetworkX directed graph."""
        G = nx.DiGraph()
        row_col = _sgc.cartesian_pair()
        for a, v in enumerate(self.vertices):
            _ = _sgc.libocn.sg_lin_to_cart(a, byref(row_col), int(self.m), int(self.n))
            pos = row_col.i, row_col.j
            G.add_node(a, pos=pos, drained_area=v["drained_area"])
        for a, v in enumerate(self.vertices):
            if v["downstream"] != _sgc.IS_ROOT:
                G.add_edge(a, v["adown"])
        return G
    
    def plot_streamgraph(self, ax=None):
        """Plot the StreamGraph in matplotlib"""
        G = self.to_networkx()
        pos = nx.get_node_attributes(G, 'pos')
        drained_area = np.asarray(list(nx.get_node_attributes(G, 'drained_area').values()))
        drained_area = (drained_area - drained_area.min()) / (drained_area.max() - drained_area.min())
        size = 50 + drained_area*1000
        lw = 1 + drained_area*5
        if ax is None:
            _, ax = plt.subplots()
        nx.draw_networkx(G, pos=pos, node_size=size, width=lw, with_labels=False, arrowstyle=ArrowStyle("-|>", head_length=1.5, head_width=0.5), ax=ax)
        return ax


def display_streamgraph(sg:StreamGraph, use_utf8=True):
    """Plot the StreamGraph in a human-readable ASCII format."""
    _sgc.libocn.sg_display(byref(sg._c_graph), use_utf8)

