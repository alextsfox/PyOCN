"""
ocn.py

High-level Python interface for working with Optimized Channel Networks (OCNs) from the libocn C library.

Author: Alexander S Fox
Copyright: (c) 2025 Alexander S Fox. All rights reserved.

This file is part of the PyOCN project.
"""

from itertools import product
import warnings

from ctypes import byref
import ctypes
from numbers import Number
from typing import Literal, Any, Callable
import networkx as nx 

import numpy as np
from tqdm import tqdm


from ._statushandler import check_status
from .utils import create_cooling_schedule, net_type_to_dag
from . import _libocn_bindings as _bindings
from . import _streamgraph_convert as sgconv
        
class OCN():
    def __init__(self, dag:nx.DiGraph, gamma:float=0.5, random_state=None, verbosity:int=0):
        """
        Main class for working with Optimized Channel Networks (OCNs). 
        Provides a high-level interface to the `libocn` C library.

        Initialize a new OCN through one of the following constructors:
            OCN.from_net_type(net_type, dims, gamma, annealing_schedule, random_state): initialize from a predefined network type.
            OCN.from_digraph(dag, gamma, annealing_schedule, random_state): initialize from an nx.DiGraph instance.

        All constructors take an initilization structure (net_type + dims, or dag) and the following:
            - gamma (float): The gamma parameter for energy calculations. Default is 0.5.
            - random_state (int, np.random.Generator, or None): Any legal input to np.random.default_rng. Used to seed the internal random number generator. Default is None (unpredictable).
            - verbosity (int): Level of verbosity for libocn output. Ranges from 0-2. Default is 0 (no output). Higher values produce more output.

        Refer to the documentation for each constructor for details on the initialization structures.
        """
        
        # validate gamma, annealing schedule, and random_state
        if not isinstance(gamma, Number):
            raise TypeError(f"gamma must be a scalar. Got {type(gamma)}.")
        if not (0 <= gamma <= 1):
            raise ValueError(f"gamma must be in the range [0, 1]. got {gamma}")
        self.gamma = gamma
        
        rng = np.random.default_rng(random_state)
        seed = rng.integers(0, int(2**32 - 1))
        _bindings.libocn.rng_seed(seed)

        # instantiate the StreamGraph_C and assign an initial energy.
        self.verbosity = verbosity
        self.__p_c_graph = sgconv.from_digraph(dag, verbose=(verbosity > 1))
        self.__p_c_graph.contents.energy = self.compute_energy()

    @classmethod
    def from_net_type(cls, net_type:str, dims:tuple, gamma=0.5, random_state=None, verbosity:int=0):
        """
        Initialize from a predefined network type and dimensions.
        
        Parameters:
            net_type (str): The type of network to create. Must be one of TODO:ALLOWED TYPES.
                Descriptions of allowed types:
                    .....
            dims (tuple): A tuple of two positive even integers specifying the dimensions of the network (rows, columns).
            gamma (float)
            random_state (int, np.random.Generator, or None)
        """
        if not isinstance(dims, tuple): raise TypeError(f"dims must be a tuple of two positive even integers, got {type(dims)}")
        if not (
            len(dims) == 2 
            and all(isinstance(d, int) and d > 0 and d % 2 == 0 for d in dims)
        ):
            raise ValueError(f"dims must be a tuple of two positive even integers, got {dims}")
        
        if verbosity == 1: print(f"Creating {net_type} network DiGraph with dimensions {dims}...", end="")
        dag = net_type_to_dag(net_type, dims)
        if verbosity == 1: print(" Done.")
        return cls(dag, gamma, random_state, verbosity=verbosity)

    @classmethod
    def from_digraph(cls, dag:nx.DiGraph, gamma=0.5, random_state=None, verbosity:int=0):
        """
        Initialize from an existing NetworkX DiGraph.
        Parameters:
            dag (nx.DiGraph): An existing directed acyclic graph (DAG) representing the stream network. Be a valid DAG for a StreamGraph (see StreamGraph for details). Additionally, must have even dimensions and at least 4 vertices.
            gamma (float)
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
        return cls(dag, gamma, random_state, verbosity=verbosity)

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
            self.gamma, 
            temperature
        ))

    def fit(
        self,
        n_iterations:int=None,
        cooling_rate:float=1.0,
        constant_phase:float=0.0,
        pbar:bool=True,
        report_energy_interval:int=1000,
    ) -> np.ndarray:
        """
        Optimize the OCN using simulated annealing. This method simulated n_iterations erosion events, using a simulated annealing algorithm to accept or reject each event.

        
        ## Simulated annealing algorithm details:

        At each iteration, a single erosion event is proposed, and the change in energy $\Delta E$ is computed.
        Erosion events are accepted or rejected with probability:

        $$
        P(accept) = e^{-\Delta E / T}
        $$

        Where $T$ is the current temperature, which is held constant for a short time after initialization, before decaying exponentially:

        $$
        T(i) =
        \begin{cases}
        E_0 & \text{if } i < C \cdot N \\
        E_0 \cdot e^{i - C \cdot N} & \mathrm{if } i \geq C \cdot N
        \end{cases}
        $$

        Where $E_0$ is the initial energy of the OCN, $N$ is the total number of iterations, and $i$ is the current iteration.

        Note that changes that decrease energy are always accepted, since $P(accept) > 1$ for $\Delta E < 0$.

        To ensure convergence, it is recommended to use a cooling rate between 0.5 and 10 and a constant phase <= 0.3. 
        A slow cooling rate and long constant phase can cause the network configuration to depart more significantly from the initial state. 
        If cooling_rate < 0.5 and constant_phase > 0.1 are used, it is suggested to increase n_iterations with respect to the default value in order to guarantee convergence.


        Parameters:
            n_iterations (int, optional): Total number of iterations to perform. Default is 40 * rows * columns.
            cooling_rate (float): The cooling rate for the annealing schedule. Must be in [0, 1]. Default is 1.0.
            constant_phase (float): The fraction of the total iterations to spend at constant temperature before cooling begins. Must be in [0, 1]. Default is 0.0 (start cooling immediately).
            pbar (bool): Whether to display a progress bar. Default is True.
            max_iterations_per_loop (int): Maximum number of iterations to perform per block. Default is 2,500. This limits memory usage when optimizing large OCNs.
            report_energy_interval (int): Interval at which to report energy in the progress bar. Default is 1,000.
        """
        max_iterations_per_loop = report_energy_interval

        if n_iterations is None:
            n_iterations = int(40*self.dims[0]*self.dims[1])
        if not (isinstance(n_iterations, int) and n_iterations > 0):
            raise ValueError(f"n_iterations must be a positive integer, got {n_iterations}")
        
        cooling_schedule = create_cooling_schedule(
            ocn=self,
            constant_phase=constant_phase,
            n_iterations=n_iterations,
            cooling_rate=cooling_rate,
        )

        if (cooling_rate < 0.5 or constant_phase > 0.1) and n_iterations <= 50*self.dims[0]*self.dims[1]:
            warnings.warn("Using cooling_rate < 0.5 and constant_phase > 0.1 with may cause convergence issues. Consider increasing n_iterations.")

        completed_iterations = 0
        energy_out = np.empty(n_iterations//report_energy_interval + 1, dtype=np.float64)
        energy_out[0] = self.energy
        
        energy_buf = np.empty(max_iterations_per_loop, dtype=np.float64)
        energy_ptr = energy_buf.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

        anneal_buf = np.empty(max_iterations_per_loop, dtype=np.float64)
        anneal_ptr = anneal_buf.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        with tqdm(total=n_iterations, desc="OCN Optimization", unit_scale=True, dynamic_ncols=True, disable=not (pbar or self.verbosity >= 1)) as pbar:
            pbar.set_postfix({"Energy": self.energy, "P(Accept)": 1.0})
            pbar.update(0)
            
            while completed_iterations < n_iterations:
                iterations_this_loop = min(max_iterations_per_loop, n_iterations - completed_iterations)
                anneal_buf[:iterations_this_loop] = cooling_schedule(
                    np.arange(completed_iterations, completed_iterations + iterations_this_loop)
                )
                check_status(_bindings.libocn.ocn_outer_ocn_loop(
                    self.__p_c_graph, 
                    energy_ptr,
                    iterations_this_loop, 
                    self.gamma, 
                    anneal_ptr
                ))
                completed_iterations += iterations_this_loop

                # copy energies to output array
                s = slice(
                    (completed_iterations - iterations_this_loop)//report_energy_interval + 1,
                    completed_iterations//report_energy_interval + 1,
                    1
                )
                energy_out[s] = energy_buf[:iterations_this_loop:report_energy_interval]

                pbar.set_postfix({"Energy": self.energy, "T": anneal_buf[iterations_this_loop - 1]})
                pbar.update(iterations_this_loop)
        return energy_out
        

__all__ = [
    "OCN",
]
