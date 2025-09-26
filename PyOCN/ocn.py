"""
ocn.py

High-level Python interface for working with Optimized Channel Networks (OCNs) from the libocn C library.

Author: Alexander S Fox
Copyright: (c) 2025 Alexander S Fox. All rights reserved.

This file is part of the PyOCN project.
"""

from itertools import product

from ctypes import byref
import ctypes
from numbers import Number
from typing import Literal, Any, Callable
import networkx as nx 

import numpy as np
from tqdm import tqdm

from ._statushandler import check_status
from . import _libocn_bindings as _bindings
from . import _streamgraph_convert as sgconv

_allowed_net_types = {"I", "H", "V", "T"}

#TODO: add ability to move root?
def net_type_to_dag(net_type:Literal["I", "H", "V", "T"], dims:tuple) -> nx.DiGraph:
    """Create a NetworkX DiGraph representing a predefined network type and dimensions.
    Parameters:
        net_type (str): The type of network to create. Must be one of "I", "H", "V", "T".
            Descriptions of allowed types:
                "I": 

                    O--O--O--O--O 
                          |
                    O--O--O--O--O
                          |
                    O--O--O--O--O
                          |
                    O--O--X--O--O

                "V": 
                    O  O  O  O  O 
                     \  \ | /  /
                    O  O  O  O  O
                     \  \ | /  /
                    O  O  O  O  O
                     \  \ | /  /
                    O--O--X--O--O

                "H": 
                    O  O  O  O 
                    |  |  | /
                    O  O  O--O
                    |  | /   
                    O  O--O--O
                    | /     
                    X--O--O--O



                "T": Channels flowing towards the center from three corners.
        dims (tuple): A tuple of two positive even integers specifying the dimensions of the network (rows, columns).
    """
    rows, cols = dims
    G = nx.DiGraph()
    match net_type:
        case "I":
            jroot = cols // 2
            for i, j in product(range(rows), range(cols)):
                n = i*cols + j
                G.add_node(n, pos=(i, j))
                if j < jroot:
                    G.add_edge(n, n+1)
                elif j > jroot:
                    G.add_edge(n, n-1)
                elif i > 0:
                    G.add_edge(n, n - cols)

        case "V":
            jroot = cols // 2
            for i, j in product(range(rows), range(cols)):
                n = i*cols + j
                G.add_node(n, pos=(i, j))
                if i > 0:
                    if j < jroot:
                        G.add_edge(n, n - cols + 1)
                    elif j > jroot:
                        G.add_edge(n, n - cols - 1)
                    else:
                        G.add_edge(n, n - cols)
                else:
                    if j < jroot:
                        G.add_edge(n, n + 1)
                    elif j > jroot:
                        G.add_edge(n, n - 1)
        case "H": # hip roof is like V, but flowing towards a corner.
            for i, j in product(range(rows), range(cols)):
                n = i*cols + j
                G.add_node(n, pos=(i, j))
                if i == j and i > 0:  # main diagonal
                    G.add_edge(n, n - cols - 1)
                elif i > j:
                    G.add_edge(n, n - cols)
                elif j > i:
                    G.add_edge(n, n - 1)
        case "T":
            raise NotImplementedError("T net type not yet implemented.")
        
    return G
        

class OCN():
    def __init__(self, dag:nx.DiGraph, gamma:float=0.5, annealing_schedule:float|Callable=None, random_state=None):
        """
        Main class for working with Optimized Channel Networks (OCNs). 
        Provides a high-level interface to the `libocn` C library.

        Initialize a new OCN through one of the following constructors:
            OCN.from_net_type(net_type, dims, gamma, annealing_schedule, random_state): initialize from a predefined network type.
            OCN.from_digraph(dag, gamma, annealing_schedule, random_state): initialize from an nx.DiGraph instance.

        All constructors take an initilization structure (net_type + dims, or dag) and the following:
            - gamma (float): The gamma parameter for energy calculations. Default is 0.5.
            - annealing_schedule (float, callable): Either (a) A function that takes the iteration number as the first positional argument and returns the temperature or (b) a constant temperature, between 0 and 1. Default is 0.0 (only accept energy-lowering moves).
            - random_state (int, np.random.Generator, or None): Any legal input to np.random.default_rng. Used to seed the internal random number generator. Default is None (unpredictable).

        Refer to the documentation for each constructor for details on the initialization structures.
        """
        
        # validate gamma, annealing schedule, and random_state
        if not isinstance(gamma, Number):
            raise TypeError(f"gamma must be a scalar. Got {type(gamma)}.")
        if not (0 <= gamma <= 1):
            raise ValueError(f"gamma must be in the range [0, 1]. got {gamma}")
        self.gamma = gamma

        self.annealing_schedule = lambda _: 0.0
        if (
            not isinstance(annealing_schedule, Number) 
            and not callable(annealing_schedule)
        ):
            raise TypeError(f"annealing_schedule must be a callable or a scalar, got {type(annealing_schedule)}")
        if isinstance(annealing_schedule, Number):
            if not (0 <= annealing_schedule <= 1):
                raise ValueError(f"If annealing_schedule is a scalar, it must be in the range [0, 1], got {annealing_schedule}")
            self.annealing_schedule = lambda _: annealing_schedule
        elif callable(annealing_schedule):
            self.annealing_schedule = annealing_schedule
        
        rng = np.random.default_rng(random_state)
        seed = rng.integers(0, int(2**32 - 1))
        _bindings.libocn.rng_seed(seed)

        # instantiate the StreamGraph_C and assign an initial energy.
        self.__p_c_graph = sgconv.from_digraph(dag)
        self.__p_c_graph.contents.energy = self.compute_energy()

    @classmethod
    def from_net_type(cls, net_type:str, dims:tuple, gamma=0.5, annealing_schedule=0.0, random_state=None):
        """
        Initialize from a predefined network type and dimensions.
        
        Parameters:
            net_type (str): The type of network to create. Must be one of TODO:ALLOWED TYPES.
                Descriptions of allowed types:
                    .....
            dims (tuple): A tuple of two positive even integers specifying the dimensions of the network (rows, columns).
            gamma (float)
            annealing_schedule (float, callable)
            random_state (int, np.random.Generator, or None)
        """
        if not isinstance(net_type, str) or net_type not in _allowed_net_types:
            raise ValueError(f"net_type must be one of {_allowed_net_types}, got {net_type}")
        if not isinstance(dims, tuple): raise TypeError(f"dims must be a tuple of two positive even integers, got {type(dims)}")
        if not (
            len(dims) == 2 
            and all(isinstance(d, int) and d > 0 and d % 2 == 0 for d in dims)
        ):
            raise ValueError(f"dims must be a tuple of two positive even integers, got {dims}")
        dag = net_type_to_dag(net_type, dims)
        return cls(dag, gamma, annealing_schedule, random_state)

    @classmethod
    def from_digraph(cls, dag:nx.DiGraph, gamma=0.5, annealing_schedule=0.0, random_state=None):
        """
        Initialize from an existing NetworkX DiGraph.
        Parameters:
            dag (nx.DiGraph): An existing directed acyclic graph (DAG) representing the stream network. Be a valid DAG for a StreamGraph (see StreamGraph for details). Additionally, must have even dimensions and at least 4 vertices.
            gamma (float)
            annealing_schedule (float, callable)
            random_state (int, np.random.Generator, or None)

        Notes:
        dag must satisfy the following properties:
            * Is a directed acyclic graph.
            * Each node has at least one attribute, `pos`, which is a non-negative (row:int, col:int) tuple giving the location of the node in a grid.
            * Must be a spanning tree over a dense grid of size (m x n). ie. each cell in the grid has exactly one node, each node has `out_degree=1` except the root node, which has `out_degree=0`.
            * Each node can only connect to one of its 8 neighboring cells (cardinal or diagonal adjacency). ie. edges cannot "skip" over rows or columns.
            * Edges cannot cross each other in the row-column plane.
            * The grid dimensions (m, n) must both be even integers.
        """
        return cls(dag, gamma, annealing_schedule, random_state)

    def __repr__(self):
        #TODO: too verbose?
        return f"<PyOCN.OCN object at 0x{id(self):x} with StreamGraph_C at 0x{ctypes.addressof(self.__p_c_graph.contents):x} and Vertex_C array at 0x{ctypes.addressof(self.__p_c_graph.contents.vertices):x}>"
    def __str__(self):
        return f"OCN(gamma={self.gamma}, energy={self.energy}, dims={self.dims}, root={self.root})"
    def __del__(self):
        try:
            _bindings.libocn.sg_destroy_safe(self.__p_c_graph)
            self.__p_c_graph = None
        except AttributeError:
            pass
    def __sizeof__(self):
        return (
            object.__sizeof__(self) +
            self.gamma.__sizeof__() +
            self.annealing_schedule.__sizeof__() +
            ctypes.sizeof(_bindings.StreamGraph_C) + ctypes.sizeof(_bindings.Vertex_C)*(self.dims[0]*self.dims[1])
        )
    
    def compute_energy(self) -> float:
        """
        Computes the energy of the current StreamGraph_C configuration.
        Returns:
            float: The computed energy.
        """
        dag = self.to_digraph()
        return max(nx.get_node_attributes(dag, 'energy').values())
    
    @property
    def energy(self) -> float:
        """
        The energy of the current OCN. Read-only.
        """
        return self.__p_c_graph.contents.energy
    
    @property
    def dims(self) -> tuple[int, int]:
        """The dimensions of the StreamGraph as (rows, columns). Read-only."""
        return (
            int(self.__p_c_graph.contents.dims.row),
            int(self.__p_c_graph.contents.dims.col)
        )
    @property
    def root(self) -> tuple[int, int]:
        """The (row, column) position of the root node in the StreamGraph grid. Read-only."""
        return (
            int(self.__p_c_graph.contents.root.row),
            int(self.__p_c_graph.contents.root.col)
        )

    def to_digraph(self) -> nx.DiGraph:
        """
        Construct and return a NetworkX DiGraph representation of the current StreamGraph.

        The returned DiGraph will have the following node attributes:
            - 'pos': (row, column) position of the node in the grid.
            - 'drained_area': The drained area of the node.
            - 'energy': The energy of each node.
        """
        dag = sgconv.to_digraph(self.__p_c_graph.contents)

        node_energies = dict()
        for node in nx.topological_sort(dag):
            da = dag.nodes[node]['drained_area']
            pred_energy = (node_energies[p] for p in dag.predecessors(node))
            node_energies[node] = da**self.gamma + sum(pred_energy)
        nx.set_node_attributes(dag, node_energies, 'energy')

        return dag
    
    def single_erosion_event(self, temperature:float):
        """
        #TODO: implement a max_tries
        Performs a single erosion event on the OCN.
        Parameters:
            temperature (float): The temperature parameter for the erosion event. Ranges from 0 (greedy) to 1 (always accept)
        """
        # StreamGraph *G, uint32_t *total_tries, double gamma, double temperature
        check_status(_bindings.libocn.ocn_single_erosion_event(
            self.__p_c_graph, 
            byref(ctypes.c_uint32(0)), 
            self.gamma, 
            temperature
        ))

    def outer_ocn_loop(self, n_iterations:int, max_iterations_per_loop:int=100):
        """
        Performs multiple erosion events with an annealing schedule. If n_iterations is large (>max_iterations_per_loop), this will break the processing into chunks to avoid excessive memory usage.
        Parameters:
            n_iterations (int): The number of erosion events to perform.
            max_iterations_per_loop (int): The maximum number of iterations to perform in a single call to the underlying C function. Default is 1000.
        """

        completed_iterations = 0
        with tqdm(total=n_iterations, desc="OCN Optimization") as pbar:
            pbar.set_postfix({"H": self.energy})
            pbar.update(0)
            
            while completed_iterations < n_iterations:
                iterations_this_loop = min(max_iterations_per_loop, n_iterations - completed_iterations)
                temperatures = [self.annealing_schedule(t) for t in range(completed_iterations, completed_iterations + iterations_this_loop)]
                temp_array = (ctypes.c_double * iterations_this_loop)(*temperatures)  # ctypes syntax for creating a C-compatible array from an iterable.
                check_status(_bindings.libocn.ocn_outer_ocn_loop(
                    self.__p_c_graph, 
                    iterations_this_loop, 
                    self.gamma, 
                    temp_array
                ))
                completed_iterations += iterations_this_loop
                
                pbar.set_postfix({"H": self.energy})
                pbar.update(iterations_this_loop)
        

__all__ = [
    "OCN",
]
