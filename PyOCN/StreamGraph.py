"""
#TODO: allow users to provide multiple DAGs that partition a space.

StreamGraph.py

High-level Python interface for the StreamGraph C library.

Author: Alexander S Fox
Copyright: (c) 2025 Alexander S Fox. All rights reserved.

This file is part of the OCN project.
"""

from ctypes import byref
from dataclasses import dataclass

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.patches import ArrowStyle

from . import _libocn_bindings as _bindings

graph_initializers = {
    "I": None
}

@dataclass(slots=True)
class Vertex:
    drained_area: int
    adown: int
    edges: int
    downstream: int
    visited: int


class StreamGraph:
    def __init__(self, shape:tuple[int, int]=None, init_structure:str|nx.DiGraph="toy"):
        """
        Initialize a StreamGraph object without any structure.

        Args:
            shape (tuple[int, int]): The dimensions of the graph as (rows, columns). If None, must provide nx.DiGraph init_structure.
            init_structure (str or nx.DiGraph): The initial structure of the graph. Options are:
                - "toy": Creates a small test graph with predefined structure. Overrides shape.
                - nx.DiGraph: A NetworkX directed graph to initialize the StreamGraph from. Overrides shape.

        Creates a directed graph of vertices arranged in a 2D grid, where each vertex has a drained area
        and flows downstream to another vertex. Initializes the underlying C StreamGraph structure and 
        prepares the vertex array for manipulation and analysis.
        """
        if init_structure == "toy":
            self._c_graph = _bindings.libocn.sg_make_test_graph()
        elif isinstance(init_structure, nx.DiGraph):
            self._c_graph = digraph_to_streamgraph(init_structure)
        else:    
            self._c_graph = _bindings.StreamGraph_C()
            
            if shape is None:
                raise TypeError("Must provide shape when init_structure is not 'toy' or a NetworkX DiGraph.")
            try:
                shape = tuple(int(x) for x in shape)
            except Exception as e:
                raise TypeError(f"shape must be a tuple of two integers. Got {shape}") from e
            if len(shape) != 2 or any(x <= 0 for x in shape):
                raise ValueError(f"shape must be a tuple of two positive integers. Got {len(shape)}")

            #TODO: implement other initializers
            raise NotImplementedError("Non-toy and non-DiGraph initializers not implemented yet.")

    def __getitem__(self, idx):
        return self.vertices[idx]
    def __repr__(self):
        return repr(self._c_graph)
    def __str__(self):
        return f"StreamGraph(m={self.m}, n={self.n}, i_root={self.i_root}, j_root={self.j_root}, energy={self.energy:.5f}, vertices=<{self.m * self.n} vertices>)"
    def __del__(self):
        _bindings.libocn.sg_destroy(byref(self._c_graph))
    
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
        vert_c = _bindings.Vertex_C()
        vertices = np.empty((self.m*self.n,), dtype=_bindings.vert_dtype)
        for a in range(self.m*self.n):
            _ = _bindings.libocn.sg_get_lin(self._c_graph, byref(vert_c), a)
            vertices[a] = (vert_c.drained_area, vert_c.adown, vert_c.edges, vert_c.downstream, vert_c.visited)
        return vertices
    
    def compute_energy(self) -> float:
        raise NotImplementedError("Energy computation not implemented yet.")
    
    # TODO: move to an OCN class/module
    # def single_erosion_event(self):
    #     """Perform a single erosion event on the StreamGraph."""
    #     status = _bindings.libocn.sg_single_erosion_event(byref(self._c_graph))
    #     if _bindings.STATUS_CODES[status] != "SUCCESS":
    #         raise RuntimeError(f"Failed to perform single erosion event: {_bindings.STATUS_CODES.get(status, 'Unknown error code')} ({status})")
        
    def to_networkx(self) -> nx.DiGraph:
        """Convert the StreamGraph to a NetworkX directed graph."""
        G = nx.DiGraph()
        row_col = _bindings.cartesian_pair()
        for a, v in enumerate(self.vertices):
            _ = _bindings.libocn.sg_lin_to_cart(a, byref(row_col), int(self.m), int(self.n))
            pos = row_col.i, row_col.j
            G.add_node(a, pos=pos, drained_area=v["drained_area"])
        for a, v in enumerate(self.vertices):
            if v["downstream"] != _bindings.IS_ROOT:
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
    _bindings.libocn.sg_display(byref(sg._c_graph), use_utf8)

def digraph_to_streamgraph(G:nx.DiGraph) -> _bindings.StreamGraph_C:
    """Convert a NetworkX directed graph to a StreamGraph. Graph must be a directed acyclic graph that represents a spanning tree of a set of nodes. Edges point downstream.
    Each node must have a position attribute `pos=(row, col)` where `row` and `col` are integers, representing an index into a 2d array of vertices.
    """

    """
    How to construct:
        * ensure that each node in the graph has only one immediate successor, and that there are no duplicate nodes.
        * find the maximum row and col value across all nodes, and create a 2d np array of vertices (as object dtypes) with those dimensions, where each vertex is set up as follows: drained_area=0, adown=0, edges=0, downstream=0, visited=-9999. Each vertex is a Vertex dataclass instance.
        * loop over the nodes in the graph in reverse topological order (ie we process successor nodes first) and for each node:
            * get its position attribute and use that to find the corresponding vertex in the 2d list.
            * compute the total number of predecessor nodes (including itself) and set the drained_area field of the vertex to that value.
            * if the node has a successor:
                get the successor's position attribute and use that to find the corresponding vertex in the 2d list. 
                Compute the linear index of the successor vertex using libocn.sg_cart_to_lin. Compute the cardinal direction (0=N, 1=NE, 2=E, 3=SE, row 0=northmost, col 0=westmost) from the current vertex to the successor vertex, and set correct bit of the edges field of the current vertex to 1 << direction: vertex.edges = vertex.edges ^ (1 << direction). Set the downstream field using the cardinal direction, set adown to be the linear index of the successor vertex.
            * if the node has no successor, set the downstream field of the current vertex to IS_ROOT (as a python value).
            * set the visited field of the current vertex to 0.
        * after processing all nodes, flatten the 2d array to a 1d array according to libocn.sg_cart_to_lin to convert from (row, col) -> (a). Loop over the 1d array and for each vertex, check that the visited field = 0. If not, raise a malformed graph error.
        * find the root node (the one with no successors) and get its position attribute.
        * create a StreamGraph_C instance with the appropriate fields set, and return it.
    """


