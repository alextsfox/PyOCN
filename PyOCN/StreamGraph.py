#TODO: allow users to provide multiple DAGs that partition a space.
"""
streamgraph.py

High-level Python interface for StreamGraph structures from the libocn C library.

Author: Alexander S Fox
Copyright: (c) 2025 Alexander S Fox. All rights reserved.

This file is part of the PyOCN project.
"""

from collections.abc import Generator
from ctypes import byref
from dataclasses import dataclass
import itertools
import warnings

import numpy as np
import networkx as nx


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
        Base class for working with StreamGraphs, which act as a framework for working with OCNs.
        This is a high-level python interface built on NetworkX that wraps the underlying C StreamGraph structure and associated functions.

        Args:
            shape (tuple[int, int]): The dimensions of the graph as (rows, columns). If None, must provide nx.DiGraph init_structure.
            init_structure (str or nx.DiGraph): The initial structure of the graph. Options are:
                - "I": TODO: IMPLEMENT
                - nx.DiGraph: A NetworkX directed graph to initialize the StreamGraph from. Must satisfy the following properties:
                    * Is a directed acyclic graph (ie. no cycles)
                    * Each node has at least one attribute, `pos`, which is a (row:int, col:int) tuple giving the location of the node in a grid.
                    * The graph must be a spanning tree over a dense grid of size (m x n). ie. each cell in the grid has exactly one node, each node has `out_degree=1` except the root node, which has `out_degree=0`.  TODO: relax this?
                    * Each node can only connect to one of its 8 neighboring cells (cardinal or diagonal adjacency). ie. edges cannot "skip" over rows or columns.
                    * Edges cannot cross each other in the row-column plane.

                    These constraints are checked when initializing from a DiGraph, throwing a ValueError if any are violated.
        """
        if isinstance(init_structure, nx.DiGraph):
            self._c_graph, self.DAG = _digraph_to_streamgraph_c(init_structure)
        else:    
            if shape is None:
                raise TypeError("Must provide shape when init_structure is not a NetworkX DiGraph.")
            try:
                shape = tuple(int(x) for x in shape)
            except Exception as e:
                raise TypeError(f"shape must be a tuple of two integers. Got {shape}.") from e
            if len(shape) != 2 or any(x < 2 for x in shape):
                raise ValueError(f"shape must be a tuple of two positive integers. Got {len(shape)}.")
            if any(x % 2 != 0 for x in shape):
                raise ValueError(f"Shape dimensions should be even. Got {shape}.")

            #TODO: implement other initializers
            raise NotImplementedError("non-DiGraph initializers not implemented yet.")

    def __len__(self) -> int:
        return self.size
    # def __getitem__(self, idx) -> Vertex|np.ndarray:
    #     """Get vertex or vertices by numpy-style indexing."""
    #     rows, cols = self.shape

    #     # Convert idx to row, col slices
    #     if isinstance(idx, tuple):  # multiple indices
    #         if len(idx) != 2:
    #             raise IndexError(f"Too many indices for StreamGraph: got {len(idx)}, expected 2")
    #         row_idx, col_idx = idx
    #     else:  # single index
    #         row_idx, col_idx = idx, slice(None)
        
    #     # Convert to actual row/col ranges
    #     def normalize_idx(index, max_val):
    #         if isinstance(index, slice):
    #             start = 0 if index.start is None else index.start  # slice(None, stop, step)
    #             stop = max_val if index.stop is None else index.stop  # slice(start, None, step)
    #             step = 1 if index.step is None else index.step  # slice(start, stop, None)
    #             if start < 0: start = max_val + start  # slice(-n, stop, step)
    #             if stop < 0: stop = max_val + stop  # slice(start, -n, step)
    #             return range(start, stop, step)
    #         elif isinstance(index, int):
    #             if index < 0: index = max_val + index  # -n
    #             return [index]
    #         else:
    #             raise IndexError(f"Unsupported index type: {index}, type {type(index)}")
        
    #     row_indices = normalize_idx(row_idx, rows)
    #     col_indices = normalize_idx(col_idx, cols)
        
    #     # Create result array with appropriate shape
    #     is_single_row = isinstance(row_idx, int)
    #     is_single_col = isinstance(col_idx, int)
        
    #     # Handle the various cases
    #     if is_single_row and is_single_col:
    #         # Single vertex case: sg[3, 4]
    #         r, c = row_indices[0], col_indices[0]
    #         v = _bindings.libocn.sg_get_cart_safe(byref(self._c_graph), _bindings.CartPair_C(row=r, col=c))
    #         return Vertex(v.drained_area, v.adown, v.edges, v.downstream, v.visited)
    #     else:
    #         # Create a numpy array of the appropriate shape
    #         result = np.empty((len(row_indices), len(col_indices)), dtype=np.object_)
            
    #         # Fill it with vertices
    #         for i, r in enumerate(row_indices):
    #             for j, c in enumerate(col_indices):
    #                 v = _bindings.libocn.sg_get_cart_safe(byref(self._c_graph), _bindings.CartPair_C(row=r, col=c))
    #                 result[i, j] = Vertex(v.drained_area, v.adown, v.edges, v.downstream, v.visited)
            
    #         # Convert to 1D array if needed
    #         if is_single_row:
    #             return result[0]
    #         elif is_single_col:
    #             return result[:, 0]
    #         else:
    #             return result
    def __iter__(self) -> Generator[Vertex]:
        """Iterate over the vertices in row-major order."""
        v_c = _bindings.Vertex_C()
        for r, c in itertools.product(range(self.shape[0]), range(self.shape[1])):
            check_status(_bindings.libocn.sg_get_cart_safe(
                byref(v_c), 
                byref(self._c_graph),
                _bindings.CartPair_C(row=r, col=c), 
            ))
            yield Vertex(v_c.drained_area, v_c.adown, v_c.edges, v_c.downstream, v_c.visited)
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


def display_streamgraph(sg:StreamGraph, use_utf8=True):
    """Plot the StreamGraph in a human-readable ASCII format."""
    _bindings.libocn.sg_display(byref(sg._c_graph), use_utf8)

def _digraph_to_streamgraph_c(G: nx.DiGraph) -> _bindings.StreamGraph_C, nx.DiGraph:
    """
    Convert a NetworkX directed graph into a StreamGraph_C. Called by the StreamGraph constructor.
    Returns:
        StreamGraph_C: The corresponding C StreamGraph structure.
        DAG: The validated input graph, node attributes updated to include the Vertex_C dataclass fields.
    """
    
    # TODO modify for the possibility of multiple input graphs
    # is a DAG
    if not isinstance(G, nx.DiGraph):
        raise TypeError(f"G must be a networkx.DiGraph, got {type(G)}")
    if not nx.is_directed_acyclic_graph(G):
        raise ValueError("Graph must be a DAG.")
    
    # TODO: modify for the possibility of multiple input graphs
    # pos attribute is valid
    pos = list(nx.get_node_attributes(G, "pos").values())
    if len(pos) != len(G.nodes):
        raise ValueError("All graph nodes must have a 'pos' attribute.")
    if any(
        not isinstance(p, (tuple, list))  # Check if position is a tuple or list
        or len(p) != 2  # Check if it has exactly two elements
        or not all(isinstance(x, (int, np.integer)) for x in p)  # Check if both elements are integers
        or any(x < 0 for x in p)  # Check if both elements are non-negative
        for p in pos
    ):
        raise ValueError("All graph node 'pos' attributes must be non-negative (row:int, col:int) tuples.")
    

    # TODO: modify for the possibility of multiple input graphs
    # spans a dense grid
    rows, cols = max(p[0] for p in pos) + 1, max(p[1] for p in pos) + 1
    for r, c in itertools.product(range(rows), range(cols)):
        if (r, c) not in pos:
            raise ValueError(f"Graph does not cover a dense {rows}x{cols} grid. Missing position ({r}, {c}).")
    if len(G.nodes) != rows * cols:
        raise ValueError(f"Graph does not cover a dense {rows}x{cols} grid (expected {rows*cols} nodes, got {len(G.nodes)}).")

    # TODO modify for the possibility of multiple input graphs
    # is a spanning tree
    if any(G.out_degree(u) > 1 for u in G.nodes):
        raise ValueError("Graph must be a spanning tree (each node has out_degree <= 1).")
    roots = [u for u in G.nodes if G.out_degree(u) == 0]
    if len(roots) != 1:
        raise ValueError(f"Graph must be a spanning tree (graph must have exactly one root (out_degree=0)). Found {len(roots)}.")
    
    # TODO modify for the possibility of multiple input graphs
    # edges only connect to adjacent nodes (no skipping)
    for u, v in G.edges:
        r1, c1 = G.nodes[u]["pos"]
        r2, c2 = G.nodes[v]["pos"]
        if max(abs(r1 - r2), abs(c1 - c2)) != 1:
            raise ValueError(f"Edge ({u}->{v}) connects non-adjacent nodes at positions {(r1,c1)} and {(r2,c2)}.")
    

    # TODO modify for the possibility of multiple input graphs
    # compute the drained area, adown, edges, downstream, and visited attributes for each node.
    # checking for crosses can come later: easier with edges defined.
    for n in nx.topological_sort(G):
        G.nodes[n]["visited"] = 0

        succs = list(G.successors(n))
        preds = list(G.predecessors(n))
        neighbors = succs + preds
        if len(neighbors) > 8 or len(neighbors) < 1:
            raise ValueError(f"Node {n} at position {G.nodes[n]['pos']} has {len(neighbors)} neighbors, but must have between 1 and 8.")
        if len(succs) > 1:
            raise ValueError(f"Node {n} at position {G.nodes[n]['pos']} has {len(succs)} successors, but must have at most 1.")

        G.nodes[n]["drained_area"] = 1 + sum(G.nodes[p]["drained_area"] for p in preds)

        def direction_bit(pos1, pos2):
            r1, c1 = pos1
            r2, c2 = pos2
            dr, dc = r2 - r1, c2 - c1
            match dr, dc:
                case -1,  0: return 0  # N
                case -1,  1: return 1  # NE
                case  0,  1: return 2  # E
                case  1,  1: return 3  # SE
                case  1,  0: return 4  # S
                case  1, -1: return 5  # SW
                case  0, -1: return 6  # W
                case -1, -1: return 7  # NW
                case _: raise ValueError(f"Nodes at positions {pos1} and {pos2} are not adjacent.")
        G.nodes[n]["edges"] = 0
        for nbr in neighbors:
            G.nodes[n]["edges"] |= (1 << direction_bit(G.nodes[n]["pos"], G.nodes[nbr]["pos"]))
        
        G.nodes[n]["downstream"] = _bindings.IS_ROOT
        G.nodes[n]["adown"] = rows * cols  # invalid index
        if len(succs):
            nsucc = succs[0]
            G.nodes[n]["downstream"] = direction_bit(G.nodes[n]["pos"], G.nodes[nsucc]["pos"])
            G.nodes[n]["adown"] = _bindings.libocn.cart_to_lin(
                _bindings.CartPair_C(row=G.nodes[nsucc]["pos"][0], col=G.nodes[nsucc]["pos"][1]),
                _bindings.CartPair_C(row=rows, col=cols)
            )

    # TODO modify for the possibility of multiple input graphs
    # check that edges do not cross each other
    for n in G.nodes:
        r, c = G.nodes[n]["pos"]
        succs = list(G.successors(n))
        if len(succs) == 0: continue  # skip root node
        match G.nodes[n]["downstream"]:
            case 1: r_check, c_check = r - 1, c      # NE flow: check N vertex
            case 7: r_check, c_check = r - 1, c      # NW flow: check N vertex
            case 3: r_check, c_check = r + 1, c      # SE flow: check S vertex
            case 5: r_check, c_check = r + 1, c      # SW flow: check S vertex
            case _: continue  # Not a diagonal flow, cannot cross

        # find the node with pos (r_check, c_check)
        cross_check_nodes = [u for u in G.nodes if G.nodes[u]["pos"] == (r_check, c_check)]
        if len(cross_check_nodes) == 0:
            raise ValueError(f"Node at position {(r_check, c_check)} does not exist!")
        if len(cross_check_nodes) > 1:
            raise ValueError(f"Multiple nodes found at position {(r_check, c_check)}!")
        cross_check_node = cross_check_nodes[0]
        cross_edges = G.nodes[cross_check_node]["edges"]
        if (
            G.nodes[n]["downstream"] == 1 and (cross_edges & (1 << 3))  # NE flow: N vertex has SE edge
            or (G.nodes[n]["downstream"] == 7 and (cross_edges & (1 << 5)))  # NW flow: N vertex has SW edge
            or (G.nodes[n]["downstream"] == 3 and (cross_edges & (1 << 1)))  # SE flow: S vertex has NE edge
            or (G.nodes[n]["downstream"] == 5 and (cross_edges & (1 << 7)))  # SW flow: S vertex has NW edge
        ):
            raise ValueError(f"Edge ({n}->{succs[0]}) crosses edge from node at position {(r_check, c_check)}.")

    # By now, the graph is validated and has all necessary attributes to create the C StreamGraph structure.
    c_graph = _bindings.StreamGraph_C()
    
    # CartPair dims = sg->dims;
    # Vertex vert = sg->vertices[a];
    # clockhand_t down_new = vert.downstream;

    # Vertex cross_check_vert;
    # CartPair check_row_col = sg_lin_to_cart(a, dims);
    
    # // Check appropriate vertex based on flow direction
    # if (down_new == 1) {
    #     check_row_col.row -= 1;  // NE flow: check N vertex
    # } else if (down_new == 7) {
    #     check_row_col.row -= 1;  // NW flow: check N vertex
    # } else if (down_new == 3) {
    #     check_row_col.row += 1;  // SE flow: check S vertex
    # } else if (down_new == 5) {
    #     check_row_col.row += 1;  // SW flow: check S vertex
    # } else {
    #     return false; // Not a diagonal flow
    # }

    # sg_get_cart_safe(&cross_check_vert, sg, check_row_col);
    
    # // Check for crossing streams
    # if (down_new == 1 && (cross_check_vert.edges & (1u << 3))) return true;  // NE flow: N vertex has SE edge
    # if (down_new == 7 && (cross_check_vert.edges & (1u << 5))) return true;  // NW flow: N vertex has SW edge
    # if (down_new == 3 && (cross_check_vert.edges & (1u << 1))) return true;  // SE flow: S vertex has NE edge
    # if (down_new == 5 && (cross_check_vert.edges & (1u << 7))) return true;  // SW flow: S vertex has NW edge
    
    # return false;


    









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
        preds = list(G.predecessors(u))
        if succs:
            vertex_down = succs[0]
            r2, c2 = positions[vertex_down]
            clockhand = get_clockhand_from_carts(r, c, r2, c2)
            vertex.downstream = clockhand
            vertex.adown = cart_to_lin(r2, c2)
        else:
            vertex.downstream = _bindings.IS_ROOT

        neighbors = succs + preds
        if neighbors:
            for v_neighbor in neighbors:
                r2, c2 = positions[v_neighbor]
                clockhand = get_clockhand_from_carts(r, c, r2, c2)
                vertex.edges |= (1 << clockhand)  # set the bit for this direction
        else:
            raise ValueError(f"Node {u} at position {(r,c)} has no neighbors; isolated nodes are not allowed.")

        
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
    # to avoid memory leaks: call sg_destroy_safe on c_graph if an exception occurs before the function returns
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            check_status(status)

        # Flatten using cart_to_lin: streamgraph may not be stored in row-major order. Using the provided function ensures consistency.
        vert_c = _bindings.Vertex_C()
        coords = _bindings.CartPair_C()
        for r, c in itertools.product(range(m), range(n)):
            a = cart_to_lin(r, c)
            v = Vmat[r][c]
            vert_c.drained_area = v.drained_area
            vert_c.adown = v.adown
            vert_c.edges = v.edges
            vert_c.downstream = v.downstream
            vert_c.visited = v.visited
            coords.row, coords.col = r, c
            status = _bindings.libocn.sg_set_cart_safe(byref(c_graph), vert_c, coords)
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
            
    except Exception as e:
        _bindings.libocn.sg_destroy_safe(byref(c_graph))
        raise e
        
    return c_graph

def validate_streamgraph(sg:StreamGraph) -> bool|str:
    """
    Validate the integrity of a StreamGraph.

    Checks:
      * Each vertex's downstream pointer is valid (either IS_ROOT or points to an adjacent vertex).
      * The graph is acyclic and has a single root.
      * The drained_area values are consistent with the graph structure.

    Parameters:
        sg (StreamGraph): The StreamGraph to validate.

    Returns:
        either True if valid, or an error message string if invalid.
    """
    try:
        G = sg.to_networkx()
        c_graph = _digraph_to_streamgraph_c(G)
        _bindings.libocn.sg_destroy_safe(byref(c_graph))
    except Exception as e:  # _digraph_to_streamgraph_c will destroy c_graph on failure
        return str(e)
    return True

__all__ = [
    "Vertex",
    "StreamGraph",
    "display_streamgraph",
    "validate_streamgraph",
]
