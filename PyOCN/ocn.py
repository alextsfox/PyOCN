"""
ocn.py

High-level Python interface for working with Optimized Channel Networks (OCNs) from the libocn C library.

Author: Alexander S Fox
Copyright: (c) 2025 Alexander S Fox. All rights reserved.

This file is part of the PyOCN project.
"""

from ctypes import byref
import ctypes
import networkx as nx

import numpy as np
from tqdm import tqdm

from ._statushandler import check_status
from . import _libocn_bindings as _bindings
from .streamgraph import StreamGraph



class OCN():
    def __init__(self, sg:StreamGraph=None, shape:tuple[int, int]=None, init_structure:str|nx.DiGraph="toy", gamma:float=0.5, annealing_schedule:callable=None, random_state=None):
        """
        Main class for working with Optimized Channel Networks (OCNs).
        Can be initialized from an existing StreamGraph, or by specifying a shape and initial structure.
        Parameters:
            sg (StreamGraph): An existing StreamGraph to initialize from. Overrides shape and init_structure if provided.
            shape (tuple[int, int]): The shape of the graph (rows, columns).
            init_structure (str | nx.DiGraph): The initial structure of the graph.
            gamma (float): The gamma parameter for energy calculations. Default is 0.5.
            annealing_schedule (callable): A function that takes the iteration number and returns the temperature. If None, a constant temperature of 0 is used.
            random_state: Seed or random number generator for reproducibility. If None, uses default randomness. Can be any legal seed argument to np.random.default_rng().

            See StreamGraph documentation for details on supported initial structures.
        """

        if sg is not None:
            if not isinstance(sg, StreamGraph):
                raise TypeError(f"sg must be a StreamGraph instance, got {type(sg)}")
            self.sg = sg
        else:
            self.sg = StreamGraph(shape=shape, init_structure=init_structure)

        self.gamma = gamma
        self.annealing_schedule = annealing_schedule if annealing_schedule is not None else lambda t: 0.0
        self.sg._c_graph.energy = self.compute_energy()
        
        rng = np.random.default_rng(random_state)
        seed = rng.integers(0, 2**32 - 1)
        _bindings.libocn.rng_seed(seed)


    def __len__(self) -> int:
        return self.sg.size
    def __repr__(self):
        return repr(self.sg._c_graph)
    def __str__(self):
        return f"StreamGraph(dims={self.sg.dims}, root={self.sg.root}, energy={self.energy}, vertices=<{self.sg.dims[0] * self.sg.dims[1]} vertices>)"
    def __del__(self):
        del self.sg

    def compute_energy(self) -> float:
        """
        Computes the energy of the current OCN configuration.
        Returns:
            float: The computed energy.
        """
        return np.sum(v.drained_area**self.gamma for v in self.sg._vertices)
    
    @property
    def energy(self) -> float:
        """
        Computes the energy of the current OCN configuration.
        Returns:
            float: The computed energy.
        """
        return self.sg._c_graph.energy
    
    
    def single_erosion_event(self, temperature:float):
        """
        Performs a single erosion event on the OCN.
        Parameters:
            temperature (float): The temperature parameter for the erosion event. Ranges from 0 (greedy) to 1 (always accept)
        """
        check_status(_bindings.libocn.ocn_single_erosion_event(byref(self.sg._c_graph), self.gamma, temperature))

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
                    byref(self.sg._c_graph), 
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
