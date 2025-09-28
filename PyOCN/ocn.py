"""
High-level Optimized Channel Network (OCN) interface.

This module provides a high-level interface to the underlying
``libocn`` C library. The :class:`OCN` class can be used for 
constructing and optimizing river network models using 
simulated annealing.

Notes
-----
- The underlying data structure managed by :class:`OCN` is a FlowGrid owned by
    ``libocn``. Pointer lifetime and destruction are handled safely within this
    class.
- Many operations convert to a NetworkX ``DiGraph`` for convenience. These
    conversions are slow, and are intended for inspection, analysis, 
    and visualization rather than tight inner loops.

See Also
--------
PyOCN.utils
        Helper functions for OCN fitting and construction
PyOCN.plotting
        Helper functions for visualization and plotting
"""
#TODO: add an __copy__ method that deepcopies the underlying FlowGrid_C

import warnings
import ctypes
from numbers import Number
import networkx as nx 

import numpy as np
from tqdm import tqdm


from ._statushandler import check_status
from .utils import create_cooling_schedule, net_type_to_dag
from . import _libocn_bindings as _bindings
from . import _flowgrid_convert as fgconv
        
class OCN:
    """
    Optimized Channel Network wrapper around ``libocn``.

    Use :meth:`OCN.from_net_type` or :meth:`OCN.from_digraph` to construct an
    instance. The class manages the C-side FlowGrid pointer and offers methods
    to compute energy, perform erosion events, and optimize via simulated
    annealing.

    Parameters
    ----------
    dag : nx.DiGraph
        Directed acyclic graph (DAG) over a dense grid that defines the initial
        stream network. See :meth:`OCN.from_digraph` for the graph requirements.
    gamma : float, default 0.5
        Exponent in the energy model, typically in the range [0, 1].
    random_state : int | numpy.random.Generator | None, optional
        Seed or generator for the internal RNG used by ``libocn``. ``None``
        yields nondeterministic seeding.
    verbosity : int, default 0
        Verbosity level for underlying library output (0–2). Higher values may
        produce diagnostic prints.

    Attributes
    ----------
    energy : float
        Current energy of the network (read-only property).
    dims : tuple[int, int]
        Grid dimensions (rows, cols) of the FlowGrid (read-only property).
    root : tuple[int, int]
        Row/column coordinates of the root node in the grid (read-only
        property).
    """
    def __init__(self, dag: nx.DiGraph, gamma: float = 0.5, random_state=None, verbosity: int = 0):
        """
        Construct an :class:`OCN` from a valid NetworkX ``DiGraph``.

        Notes
        -----
        Please use the classmethods :meth:`OCN.from_net_type` or
        :meth:`OCN.from_digraph` to instantiate an OCN.
        """
        
        # validate gamma, annealing schedule, and random_state
        if not isinstance(gamma, Number):
            raise TypeError(f"gamma must be a scalar. Got {type(gamma)}.")
        if not (0 <= gamma <= 1):
            warnings.warn(f"gamma values outside of [0, 1] may not be physically meaningful. Got {gamma}.")
        self.gamma = gamma
        
        rng = np.random.default_rng(random_state)
        seed = rng.integers(0, int(2**32 - 1))
        _bindings.libocn.rng_seed(seed)

        # instantiate the FlowGrid_C and assign an initial energy.
        self.verbosity = verbosity
        self.__p_c_graph = fgconv.from_digraph(dag, verbose=(verbosity > 1))
        self.__p_c_graph.contents.energy = self.compute_energy()

    @classmethod
    def from_net_type(cls, net_type:str, dims:tuple, gamma=0.5, random_state=None, verbosity:int=0):
        """
        Create an :class:`OCN` from a predefined network type and dimensions.

        Parameters
        ----------
        net_type : str
            Name of the network template to generate. See
            :func:`~PyOCN.utils.net_type_to_dag` for supported types.
        dims : tuple[int, int]
            Grid dimensions (rows, cols). Both must be positive even integers.
        gamma : float, default 0.5
            Exponent in the energy model.
        random_state : int | numpy.random.Generator | None, optional
            Seed or generator for RNG seeding.
        verbosity : int, default 0
            Verbosity level (0–2) for underlying library output.

        Returns
        -------
        OCN
            A newly constructed instance initialized from the specified
            template and dimensions.

        Raises
        ------
        TypeError
            If ``dims`` is not a tuple.
        ValueError
            If ``dims`` is not two positive even integers.
        """
        if not isinstance(dims, tuple):
            raise TypeError(f"dims must be a tuple of two positive even integers, got {type(dims)}")
        if not (
            len(dims) == 2 
            and all(isinstance(d, int) and d > 0 and d % 2 == 0 for d in dims)
        ):
            raise ValueError(f"dims must be a tuple of two positive even integers, got {dims}")
        
        if verbosity == 1:
            print(f"Creating {net_type} network DiGraph with dimensions {dims}...", end="")
        dag = net_type_to_dag(net_type, dims)
        if verbosity == 1:
            print(" Done.")
        return cls(dag, gamma, random_state, verbosity=verbosity)

    @classmethod
    def from_digraph(cls, dag: nx.DiGraph, gamma=0.5, random_state=None, verbosity: int = 0):
        """
        Create an :class:`OCN` from an existing NetworkX ``DiGraph``.

        Parameters
        ----------
        dag : nx.DiGraph
            Directed acyclic graph (DAG) representing the stream network.
        gamma : float, default 0.5
            Exponent in the energy model.
        random_state : int | numpy.random.Generator | None, optional
            Seed or generator for RNG seeding.
        verbosity : int, default 0
            Verbosity level (0-2) for underlying library output.

        Returns
        -------
        OCN
            A newly constructed instance encapsulating the provided graph.

        Notes
        -----
        The input graph must satisfy all of the following:

        - It is a directed acyclic graph (DAG).
        - Each node has attribute ``pos=(row:int, col:int)`` specifying its
          grid position with non-negative coordinates. Any other attributes
          are ignored.
        - It forms a spanning tree over a dense grid of shape ``(m, n)``: each
          grid cell corresponds to exactly one node; each non-root node has
          ``out_degree == 1``; the root has ``out_degree == 0``.
        - Edges connect only to one of the 8 neighbors (cardinal or diagonal),
          i.e., no jumps over rows or columns.
        - Edges do not cross in the row-column plane.
        - Both ``m`` and ``n`` are even integers, and there are at least four
          vertices.
        """
        return cls(dag, gamma, random_state, verbosity=verbosity)

    def __repr__(self):
        #TODO: too verbose?
        return f"<PyOCN.OCN object at 0x{id(self):x} with FlowGrid_C at 0x{ctypes.addressof(self.__p_c_graph.contents):x} and Vertex_C array at 0x{ctypes.addressof(self.__p_c_graph.contents.vertices):x}>"
    def __str__(self):
        return f"OCN(gamma={self.gamma}, energy={self.energy}, dims={self.dims}, root={self.root})"
    def __del__(self):
        try:
            _bindings.libocn.fg_destroy_safe(self.__p_c_graph)
            self.__p_c_graph = None
        except AttributeError:
            pass
    def __sizeof__(self):
        return (
            object.__sizeof__(self) +
            self.gamma.__sizeof__() +
            ctypes.sizeof(_bindings.FlowGrid_C) + 
            ctypes.sizeof(_bindings.Vertex_C)*(self.dims[0]*self.dims[1])
        )
    
    def compute_energy(self) -> float:
        """
        Compute the current energy of the network.

        Returns
        -------
        float
            The computed energy value.

        Notes
        -----
        This constructs a temporary ``DiGraph`` view to aggregate node
        energies from drained areas using the exponent ``gamma``.
        """
        dag = self.to_digraph()
        return max(nx.get_node_attributes(dag, 'energy').values())
    
    @property
    def energy(self) -> float:
        """
        Energy of the current OCN (read-only).

        Returns
        -------
        float
            Current energy.
        """
        return self.__p_c_graph.contents.energy
    
    @property
    def dims(self) -> tuple[int, int]:
        """
        Grid dimensions of the FlowGrid (read-only).

        Returns
        -------
        tuple[int, int]
            ``(rows, cols)``.
        """
        return (
            int(self.__p_c_graph.contents.dims.row),
            int(self.__p_c_graph.contents.dims.col)
        )
    @property
    def root(self) -> tuple[int, int]:
        """
        Root position in the grid (read-only).

        Returns
        -------
        tuple[int, int]
            ``(row, col)`` coordinates of the root node.
        """
        return (
            int(self.__p_c_graph.contents.root.row),
            int(self.__p_c_graph.contents.root.col)
        )

    def to_digraph(self) -> nx.DiGraph:
        """
        Create a NetworkX ``DiGraph`` view of the current grid.

        Returns
        -------
        nx.DiGraph
            A DAG with the following node attributes per node:

            - ``pos``: ``(row, col)`` grid position
            - ``drained_area``: drained area value
            - ``energy``: cumulative energy at the node
        """
        dag = fgconv.to_digraph(self.__p_c_graph.contents)

        node_energies = dict()
        for node in nx.topological_sort(dag):
            da = dag.nodes[node]['drained_area']
            pred_energy = (node_energies[p] for p in dag.predecessors(node))
            node_energies[node] = da**self.gamma + sum(pred_energy)
        nx.set_node_attributes(dag, node_energies, 'energy')

        return dag
    
    def single_erosion_event(self, temperature:float):
        """
        Perform a single erosion event at a given temperature.

        Parameters
        ----------
        temperature : float
            Temperature parameter governing acceptance probability. Typical
            range is a fraction of ocn.energy.

        Raises
        ------
        LibOCNError
            If the underlying C routine reports an error status.
        """
        # FlowGrid *G, uint32_t *total_tries, double gamma, double temperature
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
        energy_reports:int=1000,
    ) -> np.ndarray:
        """
        Optimize the OCN using simulated annealing.

        This performs ``n_iterations`` erosion events, accepting or rejecting
        proposals according to a temperature schedule. A proposal consists of
        changing the outflow direction of a randomly selected vertex. The
        new outflow direction is chosen uniformly from the valid neighbors.
        A proposal is valid if it maintains a well-formed graph structure.

        Parameters
        ----------
        n_iterations : int, optional
            Total number of iterations. Defaults to ``40 * rows * cols``.
            Always at least ``energy_reports * 10`` (this should only matter for
            extremely small grids, where ``rows * cols < 256``).
        cooling_rate : float, default 1.0
            Exponential decay rate parameter for the annealing schedule (in [0, 1]).
        constant_phase : float, default 0.0
            Fraction of iterations to hold temperature constant before
            exponential decay begins (in [0, 1]).
        pbar : bool, default True
            Whether to display a progress bar.
        report_energy_interval : int, default 1000
            Interval (in iterations) at which energy is sampled into the
            returned array.

        Returns
        -------
        numpy.ndarray
            Array of sampled energy values with shape
            ``(n_iterations // report_energy_interval + 1,)``. The first entry
            is the initial energy before optimization begins.

        Raises
        ------
        ValueError

        Warns
        -----
        UserWarning

        Notes
        -----
        At iteration ``i``, a proposal with energy change :math:`\Delta E` is
        accepted with probability

        .. math::

            P(\text{accept}) = e^{-\Delta E / T(i)}.

        The temperature schedule is piecewise, held constant initially and then
        decaying exponentially:

        .. math::

            T(i) = \begin{cases}
                E_0 & i < C N \\
                E_0 \cdot e^{\;i - C N} & i \ge C N
            \end{cases}

        where :math:`E_0` is the initial energy, :math:`N` is the total number
        of iterations, and :math:`C` is ``constant_phase``. Decreasing-energy
        moves (:math:`\Delta E < 0`) are always accepted.
        """
        if n_iterations is None:
            n_iterations = int(40*self.dims[0]*self.dims[1])
        if not (isinstance(n_iterations, int) and n_iterations > 0):
            raise ValueError(f"n_iterations must be a positive integer, got {n_iterations}")
        n_iterations = max(energy_reports*10, n_iterations)

        report_energy_interval = n_iterations // energy_reports

        max_iterations_per_loop = report_energy_interval

        
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
