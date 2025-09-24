#TODO: allow users to provide multiple DAGs that partition a space.
"""
streamgraph.py

High-level Python interface for StreamGraph structures from the libocn C library.

Author: Alexander S Fox
Copyright: (c) 2025 Alexander S Fox. All rights reserved.

This file is part of the PyOCN project.
"""

from ctypes import byref
import ctypes
from dataclasses import dataclass
import itertools

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.patches import ArrowStyle

from ._statushandler import check_status

from . import _libocn_bindings as _bindings

#TODO: implement all initializers
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

vert_dtype = np.dtype([
    ("drained_area", np.uint32),
    ("adown", np.uint32),
    ("edges", np.uint8),
    ("downstream", np.uint8),
    ("visited", np.uint8),
])


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
            self._c_graph = _digraph_to_streamgraph_c(init_structure)
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

    def __len__(self) -> int:
        return self.size
    def __getitem__(self, idx):
        """Get vertex or vertices by linear or 2D index. Follows numpy index notation."""
        vertices = self._vertices
        vertices_2d = np.empty(self.shape, dtype=vertices.dtype)        
        for a in range(vertices.size):
            row_col = _bindings.libocn.sg_lin_to_cart(a, self._c_graph.dims)
            row, col = row_col.row, row_col.col
            vertices_2d[row, col] = vertices[a]
        return vertices_2d[idx]
    def __repr__(self):
        return repr(self._c_graph)
    def __str__(self):
        return f"StreamGraph(dims={self.dims}, root={self.root}, vertices=<{self.dims[0] * self.dims[1]} vertices>)"
    def __del__(self):
        try:
            _bindings.libocn.sg_destroy_safe(byref(self._c_graph))
        except AttributeError:
            pass

    @property
    def dims(self) -> tuple[int, int]:
        return int(self._c_graph.dims.row), int(self._c_graph.dims.col)
    @property
    def shape(self) -> tuple[int, int]:
        return int(self._c_graph.dims.row), int(self._c_graph.dims.col)
    @property
    def size(self) -> int:
        return self.shape[0] * self.shape[1]
    @property
    def root(self) -> tuple[int, int]:
        return (int(self._c_graph.root.row), int(self._c_graph.root.col))
    
    # TODO: figure out how to make this compatible with our in-memory tiling representation (not implemented currently in streamgraph.c)
    @property
    def _vertices(self) -> np.ndarray:
        """
        Get the vertices as a structured 1d NumPy array. 
        It is not recommended to index into this array directly, as the in-memory order of the underlying libocn vertex array may change.
        Use the __getitem__ method instead 
        (ie do not do sg._vertices[i]. Instead, do sg[i] or sg[i, j].)"""
        
        # # method 1: copy to numpy array (safe, but slower)
        # buf_addr = ctypes.addressof(self._c_graph.vertices.contents)
        # buf_len = self.size * ctypes.sizeof(_bindings.Vertex_C)
        # buf = (ctypes.c_char * buf_len).from_address(buf_addr)
        # vertices = np.frombuffer(
        #     buf, 
        #     dtype=np.dtype([
        #         ("drained_area", np.uint32),
        #         ("adown", np.uint32),
        #         ("edges", np.uint8),
        #         ("downstream", np.uint8),
        #         ("visited", np.uint8),
        #     ]), 
        #     count=self.size)

        # method 2: get each vertex individually (slower, but guaranteed correct order)
        vertices = np.empty(self.size, dtype=Vertex)
        for a in range(self.size):
            vert_c = _bindings.libocn.sg_get_lin(byref(self._c_graph), a)
            vertices[a] = Vertex(vert_c.drained_area, vert_c.adown, vert_c.edges, vert_c.downstream, vert_c.visited)
        return vertices
    
    # TODO: move to an OCN class/module
    # def compute_energy(self) -> float:
    #     raise NotImplementedError("Energy computation not implemented yet.")
    
    # TODO: move to an OCN class/module
    # def single_erosion_event(self):
    #     """Perform a single erosion event on the StreamGraph."""
    #     status = _bindings.libocn.sg_single_erosion_event(byref(self._c_graph))
    #     if _bindings.STATUS_CODES[status] != "SUCCESS":
    #         raise RuntimeError(f"Failed to perform single erosion event: {_bindings.STATUS_CODES.get(status, 'Unknown error code')} ({status})")
        
    def to_networkx(self) -> nx.DiGraph:
        """Convert the StreamGraph to a NetworkX directed graph."""
        G = nx.DiGraph()
        for a, v in enumerate(self._vertices):
            row_col = _bindings.libocn.sg_lin_to_cart(a, self._c_graph.dims)
            pos = row_col.row, row_col.col
            G.add_node(a, pos=pos, drained_area=v.drained_area)
        for a, v in enumerate(self._vertices):
            if v.downstream != _bindings.IS_ROOT:
                G.add_edge(a, v.adown)
        return G

def plot_streamgraph(sg:StreamGraph, ax=None):
    """Plot the StreamGraph in matplotlib"""
    return plot_positional_digraph(sg.to_networkx(), ax=ax)

def plot_positional_digraph(G:nx.DiGraph, ax=None):
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

def _digraph_to_streamgraph_c(G: nx.DiGraph) -> _bindings.StreamGraph_C:
    """
    Convert a NetworkX directed graph into a StreamGraph_C.

    Requirements:
      * G must be a directed acyclic graph (DAG) represented as a networkx.DiGraph object.
      * Each node has exactly one pos=(row, col) attribute (ints).
      * Each node has out_degree = 1 except for the single root node, which has out_degree = 0 (i.e., tree structure).
      * All (row, col) positions form a dense (m x n) grid (spanning).  #TODO: relax this?
      * Cardinal adjacency: successor must be in one of the 8 neighboring cells.
      * No crossed edges (i.e., no two edges intersect in the row-column plane).

    Parameters:
        G (nx.DiGraph): The input directed acyclic graph.
    Returns:
        StreamGraph_C: The corresponding C StreamGraph structure.
    """
    
    if not isinstance(G, nx.DiGraph):
        raise TypeError(f"G must be a networkx.DiGraph, got {type(G)}")
    if not nx.is_directed_acyclic_graph(G):
        raise ValueError("Graph must be a DAG.")

    # Validate positions
    positions = dict()
    for u in G.nodes:
        if "pos" not in G.nodes[u]:
            raise ValueError(f"Node {u} must have a 'pos' attribute.")
        pos = G.nodes[u]["pos"]
        if isinstance(pos, np.ndarray): pos = pos.flatten().tolist()
        if not (
            isinstance(pos, (tuple, list)) 
            and len(pos) == 2 and
            all(isinstance(x, (int, np.integer)) for x in pos)
        ):
            raise ValueError(f"Node {u} pos must be a (row:int, col:int) tuple. Got {pos}")
        positions[u] = G.nodes[u]["pos"]

    rows = [pos[0] for pos in positions.values()]
    cols = [pos[1] for pos in positions.values()]
    m = max(rows) + 1
    n = max(cols) + 1

    # Check spanning dense grid
    expected_count = m * n
    if len(G.nodes) != expected_count:
        raise ValueError(f"Graph does not cover a dense {m}x{n} grid (expected {expected_count} nodes, got {len(G.nodes)}).")
    if set(positions.values()) != set((r, c) for r in range(m) for c in range(n)):
        raise ValueError("Graph positions do not cover a dense grid from (0,0) to (m-1,n-1). Missing, extra, or duplicate positions detected.")

    # Build inverse mapping: (row,col)->node
    grid_nodes = {}
    for u, (r, c) in positions.items():
        grid_nodes[(r, c)] = u

    # Degree constraints (tree flowing downstream)
    roots = [u for u in G.nodes if G.out_degree(u) == 0]
    if len(roots) != 1:
        raise ValueError(f"Graph must have exactly one root (out_degree=0). Found {len(roots)}.")
    root = roots[0]
    root_r, root_c = positions[root]

    for u in G.nodes:
        if G.out_degree(u) > 1:
            raise ValueError(f"Node {u} has out_degree {G.out_degree(u)} (>1).")

    # Compute drained_area using topological order
    drained_area = dict()
    for u in nx.topological_sort(G):
        drained_area[u] = 1 + sum(drained_area[p] for p in G.predecessors(u))

    # Prepare a temporary 2D matrix of Vertex dataclass
    Vmat = [[Vertex(0, 0, 0, 0, -9999) for _ in range(n)] for _ in range(m)]

    # helper functions/objects for next steps
    # Cardinal direction mapping (dr, dc) -> bit index
    clock_offset_map = {
        (-1,  0): 0,  # N
        (-1,  1): 1,  # NE
        ( 0,  1): 2,  # E
        ( 1,  1): 3,  # SE
        ( 1,  0): 4,  # S
        ( 1, -1): 5,  # SW
        ( 0, -1): 6,  # W
        (-1, -1): 7,  # NW
    }
    def get_clockhand_from_carts(r0, c0, r1, c1):
        # gets the direction index (0-7) from (r0,c0) to (r1,c1)
        dr = r1 - r0
        dc = c1 - c0
        key = (dr, dc)
        if key not in clock_offset_map:
            raise ValueError(f"Connected nodes at positions {(r0, c0)} are not adjacent {(r1,c1)}")
        return clock_offset_map[key]
    def cart_to_lin(i, j):
        coords = _bindings.CartPair_C(row=i, col=j)     
        dims = _bindings.CartPair_C(row=m, col=n)
        return int(_bindings.libocn.sg_cart_to_lin(coords, dims))

    # Fill vertices
    for u in reversed(list(nx.topological_sort(G))):  # reverse topo not strictly required for assignments
        r, c = positions[u]
        vertex = Vmat[r][c]
        vertex.drained_area = int(drained_area[u])
        succs = list(G.successors(u))
        if succs:
            vertex_down = succs[0]
            r2, c2 = positions[vertex_down]
            clockhand = get_clockhand_from_carts(r, c, r2, c2)
            vertex.downstream = clockhand
            vertex.edges ^= (1 << clockhand)  # flip the bit for this direction
            vertex.adown = cart_to_lin(r2, c2)
        else:
            vertex.downstream = _bindings.IS_ROOT
        vertex.visited = 0

    # Validate that all grid cells correspond to a node (possibly redundant?)
    for r, c in itertools.product(range(m), range(n)):
        if Vmat[r][c].visited != 0:
            raise ValueError(f"Malformed graph: cell {(r,c)} not populated by any node.")

    # Allocate a new C streamgraph
    c_graph = _bindings.StreamGraph_C()
    status = _bindings.libocn.sg_create_empty_safe(
        byref(c_graph),
        _bindings.CartPair_C(row=root_r, col=root_c),
        _bindings.CartPair_C(row=m, col=n),
    )
    check_status(status)

    # Flatten using cart_to_lin: streamgraph may not be stored in row-major order. Using the provided function ensures consistency.
    vert_c = _bindings.Vertex_C()
    for r, c in itertools.product(range(m), range(n)):
        a = cart_to_lin(r, c)
        v = Vmat[r][c]
        vert_c.drained_area = v.drained_area
        vert_c.adown = v.adown
        vert_c.edges = v.edges
        vert_c.downstream = v.downstream
        vert_c.visited = v.visited
        status = _bindings.libocn.sg_set_lin_safe(byref(c_graph), vert_c, a)
        check_status(status)

    # Set root indices
    c_graph.root = _bindings.CartPair_C(row=root_r, col=root_c)
    
    # Initialize energy
    c_graph.energy = 0.0

    # Final validation pass: read back vertices
    # TODO: (Skip if performance is critical)
    check_vert = _bindings.Vertex_C()
    for a in range(m * n):
        status = _bindings.libocn.sg_get_lin_safe(byref(check_vert), byref(c_graph), a)
        check_status(status)
        if check_vert.visited != 0:
            raise RuntimeError(f"Post-write validation failed at a={a}.")
        
    return c_graph

__all__ = [
    "Vertex",
    "StreamGraph",
    "plot_streamgraph",
    "display_streamgraph",
    "plot_positional_digraph",
]
