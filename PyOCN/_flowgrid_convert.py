#TODO: allow users to provide multiple DAGs that partition a space.
#TODO: allow edge-wrapping (toroidal grids).
#TODO: implement "every edge vertex is a root" option
"""
Functions for converting between NetworkX graphs and FlowGrid_C structures.
"""

from ctypes import byref, POINTER
from itertools import product
import warnings

from numpy import integer as np_integer
import networkx as nx

from ._statushandler import check_status
from . import _libocn_bindings as _bindings

def to_digraph(c_graph:_bindings.FlowGrid_C) -> nx.DiGraph:
        """Convert the FlowGrid_C to a NetworkX directed graph."""

        dag = nx.DiGraph()
        vert_c = _bindings.Vertex_C()
        pbar = product(range(c_graph.dims.row), range(c_graph.dims.col))
        for r, c in pbar:
            a = _bindings.libocn.fg_cart_to_lin(
                _bindings.CartPair_C(row=r, col=c), 
                c_graph.dims
            )
            check_status(_bindings.libocn.fg_get_lin_safe(
                byref(vert_c), 
                byref(c_graph),
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

def from_digraph(G: nx.DiGraph, resolution:float=1, verbose:bool=False, validate:bool=True) -> POINTER:
    """
    Convert a NetworkX directed graph into a FlowGrid_C. Called by the OCN constructor.

    Parameters
    ----------
    G: nx.DiGraph
        The digraph object to initialize from
    resolution: float
        The sidelength of each gridcell in meters. Default 1.
    validate: bool
        Whether to check the input graph for validity before conversion. Default True. For internal use only.
    Returns:
        p_c_graph: pointer to the created C FlowGrid structure.
    """

    if verbose:
        print("Converting DiGraph to FlowGrid_C...")
        print("--------------------------------------")

    # is a DAG
    if validate:
        if not isinstance(G, nx.DiGraph):
            raise TypeError(f"G must be a networkx.DiGraph, got {type(G)}")
        if not nx.is_directed_acyclic_graph(G):
            raise ValueError("Graph must be a DAG.")

        if verbose:
            print("\tGraph is a directed acyclic graph.")

    # pos attribute is valid
    pos_dict = nx.get_node_attributes(G, "pos")
    pos = list(pos_dict.values())
    pos_set = set(pos_dict.values())
    pos_dict_reversed = {v:k for k,v in pos_dict.items()}

    if validate:
        if len(pos) != len(G.nodes):
            raise ValueError("All graph nodes must have a 'pos' attribute.")
        if any(
            not isinstance(p, (tuple, list))  # Check if position is a tuple or list
            or len(p) != 2  # Check if it has exactly two elements
            or not all(isinstance(x, (int, np_integer)) for x in p)  # Check if both elements are integers
            or any(x < 0 for x in p)  # Check if both elements are non-negative
            for p in pos
        ):
            raise ValueError("All graph node 'pos' attributes must be non-negative (row:int, col:int) tuples.")
        
        if verbose:
            print("\tGraph 'pos' attributes are valid.")

    # spans a dense grid
    rows, cols = max(p[0] for p in pos) + 1, max(p[1] for p in pos) + 1

    if validate:
        if set(product(range(rows), range(cols))).difference(pos_set):
            raise ValueError(f"Graph does not cover a dense {rows}x{cols} grid.")
        if len(G.nodes) != rows * cols:
            raise ValueError(f"Graph does not cover a dense {rows}x{cols} grid (expected {rows*cols} nodes, got {len(G.nodes)}).")

        if verbose:
            print(f"\tGraph covers a dense {rows}x{cols} grid.")

        # is a spanning tree
        if any(G.out_degree(u) > 1 for u in G.nodes):
            raise ValueError("Graph must be a spanning tree (each node has out_degree <= 1).")
        roots = [u for u in G.nodes if G.out_degree(u) == 0]
        
        if verbose:
            print(f"\tFound {len(roots)} spanning trees.")

        # edges only connect to adjacent nodes (no skipping)
        for u, v in G.edges:
            r1, c1 = G.nodes[u]["pos"]
            r2, c2 = G.nodes[v]["pos"]
            if max(abs(r1 - r2), abs(c1 - c2)) != 1:
                raise ValueError(f"Edge ({u}->{v}) connects non-adjacent nodes at positions {(r1,c1)} and {(r2,c2)}.")
        
        if verbose:
            print("\tEdges connect only to adjacent nodes.")


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
    # compute the drained area, adown, edges, downstream, and visited attributes for each node.
    # checking for crosses can come later: easier with edges defined.
    for n in nx.topological_sort(G):
        G.nodes[n]["visited"] = 0

        succs = list(G.successors(n))
        preds = list(G.predecessors(n))
        neighbors = succs + preds
        if len(neighbors) > 8:
            raise ValueError(f"Node {n} at position {G.nodes[n]['pos']} has {len(neighbors)} neighbors, but must have between 0 and 8.")
        if len(succs) > 1:
            raise ValueError(f"Node {n} at position {G.nodes[n]['pos']} has {len(succs)} successors, but must have at most 1.")

        G.nodes[n]["drained_area"] = resolution**2 + sum(G.nodes[p]["drained_area"] for p in preds)
        G.nodes[n]["edges"] = 0
        for nbr in neighbors:
            G.nodes[n]["edges"] |= (1 << direction_bit(G.nodes[n]["pos"], G.nodes[nbr]["pos"]))
        
        G.nodes[n]["downstream"] = _bindings.IS_ROOT
        G.nodes[n]["adown"] = rows * cols  # invalid index
        if len(succs):
            nsucc = succs[0]
            G.nodes[n]["downstream"] = direction_bit(G.nodes[n]["pos"], G.nodes[nsucc]["pos"])
            G.nodes[n]["adown"] = _bindings.libocn.fg_cart_to_lin(
                _bindings.CartPair_C(row=G.nodes[nsucc]["pos"][0], col=G.nodes[nsucc]["pos"][1]),
                _bindings.CartPair_C(row=rows, col=cols)
            )

    if verbose:
        print("\tComputed node attributes (drained_area, adown, edges, downstream, visited).")

    # check that edges do not cross each other
    #TODO: this is the current bottleneck of the conversion process. Consider profiling.
    if validate:
        for n in G.nodes:
            r, c = G.nodes[n]["pos"]
            
            down = G.nodes[n]["downstream"]
            if down % 2 == 0: continue  # Not a diagonal flow, cannot cross
            
            succs = list(G.successors(n))
            if len(succs) == 0: continue  # skip root node

            match down:
                case 1: r_check, c_check = r - 1, c      # NE flow: check N vertex
                case 7: r_check, c_check = r - 1, c      # NW flow: check N vertex
                case 3: r_check, c_check = r + 1, c      # SE flow: check S vertex
                case 5: r_check, c_check = r + 1, c      # SW flow: check S vertex
                case _: continue  # should have already been caught

            # find the node with pos (r_check, c_check)
            cross_check_node = pos_dict_reversed.get((r_check, c_check))
            if cross_check_node is None:
                raise ValueError(f"Node at position {(r_check, c_check)} does not exist!")
            cross_edges = G.nodes[cross_check_node]["edges"]
            if (
                G.nodes[n]["downstream"] == 1 and (cross_edges & (1 << 3))  # NE flow: N vertex has SE edge
                or (G.nodes[n]["downstream"] == 7 and (cross_edges & (1 << 5)))  # NW flow: N vertex has SW edge
                or (G.nodes[n]["downstream"] == 3 and (cross_edges & (1 << 1)))  # SE flow: S vertex has NE edge
                or (G.nodes[n]["downstream"] == 5 and (cross_edges & (1 << 7)))  # SW flow: S vertex has NW edge
            ):
                raise ValueError(f"Edge ({n}->{succs[0]}) crosses edge from node at position {(r_check, c_check)}.")

        if verbose:
            print("\tChecked for crossing edges.")

    # By now, the graph is validated and has all necessary attributes to create the C FlowGrid structure.
    p_c_graph = _bindings.libocn.fg_create_empty_safe(_bindings.CartPair_C(row=rows, col=cols))
    if not p_c_graph:
        raise MemoryError("Failed to allocate memory for FlowGrid_C.")
    for n in G.nodes:
        r, c = G.nodes[n]["pos"]
        a = _bindings.libocn.fg_cart_to_lin(
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
            check_status(_bindings.libocn.fg_set_lin_safe(p_c_graph, v_c, a))
        except Exception as e:
            _bindings.libocn.fg_destroy_safe(p_c_graph)
            p_c_graph = None
            raise e
    
    p_c_graph.contents.resolution = float(resolution)

    p_c_graph.contents.nroots = len([n for n in G.nodes if G.out_degree(n) == 0])

    if p_c_graph.contents.nroots > 1:
        warnings.warn(f"FlowGrid has {p_c_graph.contents.nroots} root nodes (nodes with no downstream). This will slow down certain operations.")
    
    # do not set energy

    if verbose:
        print("\tSuccessfully created FlowGrid_C from DiGraph.")

    return p_c_graph

def validate_digraph(dag:nx.DiGraph, verbose:bool=False) -> bool|str:
    """
    Validate the integrity of a FlowGrid.

    Parameters:
        dag (nx.DiGraph): The directed acyclic graph to validate.

    Returns:
        either True if valid, or an error message string if invalid.
    """
    try:
        p_c_graph = from_digraph(dag, verbose=verbose)
        _bindings.libocn.fg_destroy_safe(p_c_graph)
        p_c_graph = None
    except Exception as e:  # _digraph_to_flowgrid_c will destroy p_c_graph on failure
        return str(e)
    return "Graph is valid."

def validate_flowgrid(c_graph:_bindings.FlowGrid_C, verbose:bool=False) -> bool|str:
    """
    Validate the integrity of a FlowGrid_C structure.

    Parameters:
        c_graph (FlowGrid_C): The FlowGrid_C structure to validate.

    Returns:
        either True if valid, or an error message string if invalid.
    """
    try:
        dag = to_digraph(c_graph, verbose=verbose)
        p_c_graph = from_digraph(dag, verbose=verbose)
        _bindings.libocn.fg_destroy_safe(p_c_graph)
        p_c_graph = None
    except Exception as e:  # _digraph_to_flowgrid_c will destroy p_c_graph on failure
        return str(e)
    return "Graph is valid."
