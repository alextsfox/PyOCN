#TODO: allow users to provide multiple DAGs that partition a space.
"""
streamgraph.py

High-level Python interface for StreamGraph structures from the libocn C library.

Author: Alexander S Fox
Copyright: (c) 2025 Alexander S Fox. All rights reserved.

This file is part of the PyOCN project.
"""

from collections.abc import Generator
from ctypes import byref, POINTER
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

        Attributes:
            dag (nx.DiGraph): The stream graph represented as a NetworkX DiGraph. Nodes in dag have the following attributes:
                - pos (tuple[int, int]): The (row, column) position of the node in the grid.
                - drained_area (int): The total drained area (number of nodes, including itself) upstream of the node.

            root (tuple[int, int]): The (row, column) position of the root node.
            shape, dims (tuple[int, int]): The dimensions of the graph as (rows, columns).
            size (int): The total number of vertices in the graph.

        Warning: the _p_c_graph attribute is for internal use only. 
        Directly modifying the _p_c_graph attribute could cause memory corruption.
        Make updates to the StreamGraph through the .dag attribute, or by instantiating a new StreamGraph object.
        """
        if isinstance(init_structure, nx.DiGraph):
            self._p_c_graph = _digraph_to_streamgraph_c(init_structure)
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
    
    def __repr__(self):
        return repr(self._p_c_graph)
    def __str__(self):
        return f"StreamGraph(dims={self.dims}, root={self.root}, vertices=<{self.dims[0] * self.dims[1]} vertices>)"
    def __del__(self):
        try:
            _bindings.libocn.sg_destroy_safe(self._p_c_graph)
            self._p_c_graph = None
        except AttributeError:
            pass

    @property
    def dims(self) -> tuple[int, int]:
        return int(self._p_c_graph.contents.dims.row), int(self._p_c_graph.contents.dims.col)
    @property
    def shape(self) -> tuple[int, int]:
        return int(self._p_c_graph.contents.dims.row), int(self._p_c_graph.contents.dims.col)
    @property
    def size(self) -> int:
        return self.shape[0] * self.shape[1]
    @property
    def root(self) -> tuple[int, int]:
        return (int(self._p_c_graph.contents.root.row), int(self._p_c_graph.contents.root.col))
    @property
    def dag(self) -> nx.DiGraph:
        """
        The stream graph represented as a NetworkX DiGraph. 
        Internal changes to the attributes of .dag will not be reflected in the fundamental StreamGraph structure.
        To update the internal state of the StreamGraph, you should explicitly set the .dag attribute.

        If modification of the .dag throws an exception, the underlying structure will remain unchanged.

        Example:
            sg = StreamGraph(dag)
            
            # WRONG: this will *not* work as expected
            # sg.dag.add_edge(5, 10)
            
            # CORRECT:
            nx_dag = sg.dag
            nx_dag.add_edge(5, 10)
            sg.dag = nx_dag
        """
        return self._to_networkx()
    @dag.setter
    def dag(self, dag: nx.DiGraph):
        new_p_c_graph = _digraph_to_streamgraph_c(dag)
        if self._p_c_graph:
            _bindings.libocn.sg_destroy_safe(self._p_c_graph)
            self._p_c_graph = None
        self._p_c_graph = new_p_c_graph

    def _to_networkx(self) -> nx.DiGraph:
        """Convert the StreamGraph to a NetworkX directed graph."""
        dag = nx.DiGraph()
        vert_c = _bindings.Vertex_C()
        for r, c in itertools.product(range(self.shape[0]), range(self.shape[1])):
            a = _bindings.libocn.sg_cart_to_lin(
                _bindings.CartPair_C(row=r, col=c), 
                self._p_c_graph.contents.dims
            )
            check_status(_bindings.libocn.sg_get_lin_safe(
                byref(vert_c), 
                self._p_c_graph,
                a,
            ))
            dag.add_node(
                a, 
                pos=(r, c), 
                drained_area=vert_c.drained_area,
                _adown=vert_c.adown,
                _edges=vert_c.edges,
                _downstream=vert_c.downstream,
                _visited=vert_c.visited,
            )
            if vert_c.downstream != _bindings.IS_ROOT:
                dag.add_edge(a, vert_c.adown)
        
        return dag


def display_streamgraph(sg:StreamGraph, ascii=True):
    """Print the StreamGraph in a human-readable format
    
    Parameters
    ----------
    sg : StreamGraph
        The StreamGraph to display.
    ascii : bool, optional
        Whether to use ASCII characters for display (default is True). If False, uses UTF-8 (prettier).
    """
    _bindings.libocn.sg_display(sg._p_c_graph, not ascii)

def _digraph_to_streamgraph_c(G: nx.DiGraph) -> POINTER:
    """
    Convert a NetworkX directed graph into a StreamGraph_C. Called by the StreamGraph constructor.
    Returns:
        p_c_graph: The created C StreamGraph structure.
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
            G.nodes[n]["adown"] = _bindings.libocn.sg_cart_to_lin(
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
    p_c_graph = _bindings.libocn.sg_create_empty_safe(_bindings.CartPair_C(row=rows, col=cols))
    if not p_c_graph:
        raise MemoryError("Failed to allocate memory for StreamGraph_C.")
    for n in G.nodes:
        r, c = G.nodes[n]["pos"]
        a = _bindings.libocn.sg_cart_to_lin(
            _bindings.CartPair_C(row=r, col=c),
            _bindings.CartPair_C(row=rows, col=cols)
        )
        v_c = _bindings.Vertex_C(
            drained_area=G.nodes[n]["drained_area"],
            adown=G.nodes[n]["adown"],
            edges=G.nodes[n]["edges"],
            downstream=G.nodes[n]["downstream"],
            visited=G.nodes[n]["visited"],
        )
        try:
            check_status(_bindings.libocn.sg_set_lin_safe(p_c_graph, v_c, a))
        except Exception as e:
            _bindings.libocn.sg_destroy_safe(p_c_graph)
            p_c_graph = None
            raise e
        
    # set root
    p_c_graph.contents.root = _bindings.CartPair_C(row=G.nodes[roots[0]]["pos"][0], col=G.nodes[roots[0]]["pos"][1])
    
    # do not set energy
    
    return p_c_graph

def validate_streamgraph(sg:StreamGraph) -> bool|str:
    """
    Validate the integrity of a StreamGraph.

    Parameters:
        sg (StreamGraph): The StreamGraph to validate.

    Returns:
        either True if valid, or an error message string if invalid.
    """
    try:
        G = sg.to_networkx()
        p_c_graph = _digraph_to_streamgraph_c(G)
        _bindings.libocn.sg_destroy_safe(p_c_graph)
        p_c_graph = None
    except Exception as e:  # _digraph_to_streamgraph_c will destroy p_c_graph on failure
        return str(e)
    return "sg is valid."

__all__ = [
    "StreamGraph",
    "display_streamgraph",
    "validate_streamgraph",
]
