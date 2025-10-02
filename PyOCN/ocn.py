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
from typing import Any
from os import PathLike
from numbers import Number
from pathlib import Path
from dataclasses import dataclass

import networkx as nx 
import numpy as np
from tqdm import tqdm

from ._statushandler import check_status
from .utils import create_cooling_schedule, net_type_to_dag
from . import _libocn_bindings as _bindings
from . import _flowgrid_convert as fgconv

@dataclass(slots=True, frozen=True)
class FitResult:
    """
    Result of OCN optimization.
    
    Attributes
    -----------
    i_raster : np.ndarray
        List of iteration indices at which array snapshots were recorded.
        Shape (n_reports,).
    energy_rasters : np.ndarray
        Array snapshots of energy recorded during fitting, shaped
        (n_reports, band, rows, cols).
    area_rasters: np.ndarray
        Array snapshots of drained area recorded during fitting, shaped
        (n_reports, rows, cols).
    i_energy : np.ndarray
        List of iteration indices at which total energy snapshots were recorded.
        Shape (n_reports,).
    energies : np.ndarray
        Total energy snapshots recorded during fitting.
        Shape (n_reports,).
    temperatures : np.ndarray
        Temperature values at each recorded energy snapshot.
        Shape (n_reports,).
    """
    i_raster: np.ndarray
    energy_rasters: np.ndarray
    area_rasters: np.ndarray
    i_energy: np.ndarray
    energies: np.ndarray
    temperatures: np.ndarray

    def __str__(self):
        return (f"FitResult(energies={len(self.energies)}, grids={len(self.energy_rasters)})")

class OCN:
    """
    The main class for interacting with Optimized Channel Networks. 

    Use :meth:`OCN.from_net_type` or :meth:`OCN.from_digraph` to construct an
    instance. 
    
    Constructor Methods
    -------------------
    :meth:`from_net_type`
        Create an OCN from a predefined network type and dimensions.
    :meth:`from_digraph`
        Create an OCN from an existing NetworkX DiGraph.

    Export Methods
    --------------
    :meth:`to_digraph`
        Export the current grid to a NetworkX DiGraph.
    :meth:`to_numpy`
        Export raster arrays (energy, drained area) as numpy arrays.
    :meth:`to_xarray`
        Export raster arrays as an xarray Dataset (requires xarray).
    :meth:`to_gtiff`
        Export raster arrays to a GeoTIFF file (requires rasterio).

    Optimization Methods
    --------------------
    :meth:`compute_energy`
        Compute the current energy of the network.
    :meth:`single_erosion_event`
        Perform a single erosion event at a given temperature.
    :meth:`fit`
        Optimize the network using simulated annealing.

    Attributes
    ----------
    energy : float
        Current energy of the network (read-only property).
    dims : tuple[int, int]
        Grid dimensions (rows, cols) of the FlowGrid (read-only property).
    resolution: float
        The side length of each grid cell in meters (read-only property).
    nroots : int
        Number of root nodes in the current OCN grid (read-only property).
    gamma : float
        Exponent in the energy model.
    master_seed : int
        The seed used to initialize the internal RNG.
    verbosity : int
        Verbosity level for underlying library output (0-2).
    
    
    Examples
    --------
    The following is a simple example of creating, optimizing, and plotting
    an OCN using PyOCN and Matplotlib. More examples are available in the
    `demo.ipynb` notebook in the repository (https://github.com/alextsfox/PyOCN).

    >>> import PyOCN as po
    >>> import matplotlib as mpl
    >>> import matplotlib.pyplot as plt
    >>> # initialize a V-shaped network on a 64x64 grid
    >>> ocn = po.OCN.from_net_type(
    >>>     net_type="V",
    >>>     dims=(64, 64),
    >>>     random_state=8472,
    >>> )
    >>> # optimize the network with simulated annealing
    >>> result = ocn.fit(pbar=True)
    >>> # plot the energy and drained area evolution, along with the final raster
    >>> fig, axs = plt.subplots(2, 1, height_ratios=[1, 4], figsize=(6, 8))
    >>> axs[0].plot(result.i_energy, result.energies / np.max(result.energies))
    >>> axs[0].plot(result.i_energy, result.temperatures / np.max(result.temperatures))
    >>> axs[0].set_ylabel("Normalized Energy")
    >>> axs[0].set_xlabel("Iteration")
    >>> norm = mpl.colors.PowerNorm(gamma=0.5)
    >>> po.plot_ocn_raster(ocn, attribute='drained_area', norm=norm, cmap="Blues", ax=axs[1])
    >>> axs[1].set_axis_off()
    """
    def __init__(self, dag: nx.DiGraph, resolution: float=1.0, gamma: float = 0.5, random_state=None, verbosity: int = 0, validate:bool=True):
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

        self.master_seed = random_state
        

        if not isinstance(resolution, Number):
            raise TypeError(f"resolution must be numeric. Got {type(resolution)}")
        # instantiate the FlowGrid_C and assign an initial energy.
        self.verbosity = verbosity
        self.__p_c_graph = fgconv.from_digraph(dag, resolution, verbose=(verbosity > 1), validate=validate)
        self.__p_c_graph.contents.energy = self.compute_energy()

    @classmethod
    def from_net_type(cls, net_type:str, dims:tuple, resolution:float=1, gamma=0.5, random_state=None, verbosity:int=0):
        """
        Create an :class:`OCN` from a predefined network type and dimensions.

        Parameters
        ----------
        net_type : str
            Name of the network template to generate. See
            :func:`~PyOCN.utils.net_type_to_dag` for supported types.
        dims : tuple[int, int]
            Grid dimensions (rows, cols). Both must be positive even integers.
        resolution : int, optional
            The side length of each grid cell in meters.
        gamma : float, default 0.5
            Exponent in the energy model.
        random_state : int | numpy.random.Generator | None, optional
            Seed or generator for RNG seeding.
        verbosity : int, default 0
            Verbosity level (0-2) for underlying library output.

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
        return cls(dag, resolution, gamma, random_state, verbosity=verbosity, validate=False)

    @classmethod
    def from_digraph(cls, dag: nx.DiGraph, resolution:float=1, gamma=0.5, random_state=None, verbosity: int = 0):
        """
        Create an :class:`OCN` from an existing NetworkX ``DiGraph``.

        Parameters
        ----------
        dag : nx.DiGraph
            Directed acyclic graph (DAG) representing the stream network.
        resolution : int, optional
            The side length of each grid cell in meters.
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
        - The graph can be partitioned into one or more spanning trees over 
          a dense grid of shape ``(m, n)``: each grid cell corresponds to 
          exactly one node; each non-root node has ``out_degree == 1``; 
          the roots have ``out_degree == 0``.
        - Edges connect only to one of the 8 neighbors (cardinal or diagonal),
          i.e., no jumps over rows or columns.
        - Edges do not cross in the row-column plane.
        - Both ``m`` and ``n`` are even integers, and there are at least four
          vertices.
        """
        return cls(dag, resolution, gamma, random_state, verbosity=verbosity, validate=True)

    def __repr__(self):
        #TODO: too verbose?
        return f"<PyOCN.OCN object at 0x{id(self):x} with FlowGrid_C at 0x{ctypes.addressof(self.__p_c_graph.contents):x} and Vertex_C array at 0x{ctypes.addressof(self.__p_c_graph.contents.vertices):x}>"
    def __str__(self):
        return f"OCN(gamma={self.gamma}, energy={self.energy}, dims={self.dims}, resolution={self.resolution}m, verbosity={self.verbosity})"
    def __del__(self):
        try:
            _bindings.libocn.fg_destroy_safe(self.__p_c_graph)
            self.__p_c_graph = None
        except AttributeError:
            pass
    def __sizeof__(self) ->int:
        return (
            object.__sizeof__(self) +
            self.gamma.__sizeof__() +
            ctypes.sizeof(_bindings.FlowGrid_C) + 
            ctypes.sizeof(_bindings.Vertex_C)*(self.dims[0]*self.dims[1])
        )

    def __copy__(self) -> "OCN":
        """
        Create a deep copy of the OCN, including the underlying FlowGrid_C.
        Also copies the current RNG state. The new copy and the original
        will be independent from each other and behave identically statistically.

        If you want the copy to have a different random state, call :meth:`reseed`
        after copying.
        """
        cpy = object.__new__(type(self))
        cpy.gamma = self.gamma
        cpy.verbosity = self.verbosity
        cpy.master_seed = self.master_seed
        cpy.resolution = self.resolution

        cpy_p_c_graph = _bindings.libocn.fg_copy_safe(self.__p_c_graph)
        if not cpy_p_c_graph:
            raise MemoryError("Failed to copy FlowGrid_C in OCN.__copy__")
        cpy.__p_c_graph = cpy_p_c_graph
        return cpy

    def __deepcopy__(self, memo) -> "OCN":
        """
        Create a deep copy of the OCN, including the underlying FlowGrid_C.
        Also copies the current RNG state. The new copy and the original
        will be independent from each other and behave identically statistically.
        """
        return self.__copy__()

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
        return _bindings.libocn.ocn_compute_energy(self.__p_c_graph, self.gamma)
    
    @property
    def energy(self) -> float:
        """
        Energy of the current OCN.

        Returns
        -------
        float
            Current energy.
        """
        return self.__p_c_graph.contents.energy
    @property
    def resolution(self) -> float:
        """
        Resolution of the current OCN grid in m (read-only).

        Returns
        -------
        float
            Current resolution.
        """
        return self.__p_c_graph.contents.resolution
    @property
    def nroots(self) -> int:
        """
        Number of root nodes in the current OCN grid (read-only).

        Returns
        -------
        int
            Current number of root nodes.
        """
        return int(self.__p_c_graph.contents.nroots)
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
    def master_seed(self) -> int:
        """
        The seed used to initialize the internal RNG.

        Returns
        -------
        int
            Current master seed.
        """
        return self.__master_seed
    @master_seed.setter
    def master_seed(self, random_state:int|None|np.random.Generator=None):
        """
        Seed the internal RNG.

        Parameters
        ----------
        random_state : int | numpy.random.Generator | None, optional
            Seed or generator for RNG. If None, system entropy is used.
            If an integer or Generator is provided, the
            new seed is drawn from it directly.
        """
        if not isinstance(random_state, (int, type(None), np.random.Generator)):
            raise ValueError("RNG must be initialized with an integer/Generator/None.")
        self.__master_seed = np.random.default_rng(random_state).integers(0, int(2**32 - 1))
        _bindings.libocn.rng_seed(self.master_seed)
    def _advance_seed(self):
        random_state = self.master_seed
        new_random_state = np.random.default_rng(random_state).integers(0, int(2**32 - 1))
        self.__master_seed = new_random_state
        _bindings.libocn.rng_seed(self.__master_seed)

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
            node_energies[node] = (
                dag.nodes[node]['drained_area']**self.gamma 
                + sum(node_energies[p] for p in dag.predecessors(node))
            )
        nx.set_node_attributes(dag, node_energies, 'energy')

        return dag
    
    def to_gtiff(self, west:float, north:float, crs: Any, path:str|PathLike):
        """
        Export a raster of the current FlowGrid to a .gtiff file
        using rasterio. The resulting raster has 2 bands: `energy`
        and `drained_area`

        N.B. This uses the .resolution attribute to set pixel size in the raster
        This is assumed to be in units as meters. Using a CRS with different units
        may have unexpected results! It is recommended to choose a CRS with meter units
        then transform to another CRS later if needed.

        Parameters
        ----------
        west : float
            The western border of the raster in crs units, corresponding
            to column 0.
        north : float
            The western border of the raster in crs units, corresponding
            to row 0.
        crs : Any
            The crs for the resulting gtiff, passed to `rasterio.open`
        path : str or Pathlike
            The output path for the resulting gtiff file.



        Returns
        -------
        nx.DiGraph
            A DAG with the following node attributes per node:

            - ``pos``: ``(row, col)`` grid position
            - ``drained_area``: drained area value
            - ``energy``: cumulative energy at the node
        """
        try:
            import rasterio
            from rasterio.transform import from_origin
        except ImportError as e:
            raise ImportError(
                "PyOCN.OCN.to_gtiff() requires rasterio to be installed. Install with `pip install rasterio`."
            ) from e

        dag = self.to_digraph()
        energy = np.zeros(self.dims)   
        drained_area = np.zeros(self.dims)   
        for node in dag.nodes:
            r, c = dag.nodes[node]['pos']
            energy[r, c] = dag.nodes[node]['energy']
            drained_area[r, c] = dag.nodes[node]['drained_area']
        
        transform = from_origin(west, north, self.resolution, self.resolution)

        # Write two bands: 1=energy, 2=drained_area
        with rasterio.open(
            Path(path),
            "w",
            driver="GTiff",
            height=self.dims[0],
            width=self.dims[1],
            count=2,
            dtype=str(energy.dtype),
            crs=crs,
            transform=transform,
            compress="deflate",
        ) as dst:
            dst.write(energy, 1)
            dst.write(drained_area, 2)
            # Band descriptions (nice to have)
            try:
                dst.set_band_description(1, "energy")
                dst.set_band_description(2, "drained_area")
            except Exception:
                pass
    
    def to_numpy(self):
        """
        Export the current FlowGrid to a numpy array with shape (2, rows, cols).
        Has two channels: 0=energy, 1=drained_area.
        """
        dag = self.to_digraph()
        energy = np.zeros(self.dims)
        for node in dag.nodes:
            r, c = dag.nodes[node]['pos']
            energy[r, c] = dag.nodes[node]['energy']
        drained_area = np.zeros(self.dims)
        for node in dag.nodes:
            r, c = dag.nodes[node]['pos']
            drained_area[r, c] = dag.nodes[node]['drained_area']
        return np.stack([energy, drained_area], axis=0)
    
    def to_xarray(self):
        try:
            import xarray as xr
        except ImportError as e:
            raise ImportError(
                "PyOCN.OCN.to_xarray() requires xarray to be installed. Install with `pip install xarray`."
            ) from e
        array_out = self.to_numpy()
        return xr.Dataset(
            data_vars={
                "energy_rasters": (["y", "x"], array_out[0]),
                "area_rasters": (["y", "x"], array_out[1]),
            },
            coords={
                "y": ("y", np.arange(self.dims[0])*self.resolution),
                "x": ("x", np.arange(self.dims[1])*self.resolution),
            },
            attrs={
                "description": "OCN fit result arrays",
                "resolution_m": self.resolution,
                "gamma": self.gamma,
                "master_seed": int(self.master_seed),
            }
        )

    def single_erosion_event(self, temperature:float, xarray_out:bool=False) -> "FitResult | xr.Dataset":
        """
        Perform a single erosion event at a given temperature.

        Parameters
        ----------
        temperature : float
            Temperature parameter governing acceptance probability. Typical
            range is a fraction of ocn.energy.
        xarray_out : bool, default False
            If True, the returned result will be an xarray.Dataset instead of a FitResult.
            Requires xarray to be installed.

        Raises
        ------
        LibOCNError
            If the underlying C routine reports an error status.
        """
        # FlowGrid *G, uint32_t *total_tries, double gamma, double temperature
        self._advance_seed()
        check_status(_bindings.libocn.ocn_single_erosion_event(
            self.__p_c_graph,
            self.gamma, 
            temperature
        ))

        array_out = self.to_numpy()
        if not xarray_out:
            return FitResult(
                i_raster = np.array([0], dtype=np.int32),
                energy_rasters = array_out[np.newaxis, 0],
                area_rasters = array_out[np.newaxis, 1],
                i_energy = np.array([0], dtype=np.int32),
                energies = np.array([self.energy], dtype=np.float64),
                temperatures = np.array([temperature], dtype=np.float64),
            )
        try:
            import xarray as xr
        except ImportError as e:
            raise ImportError(
                "PyOCN.OCN.fit() with xarray_out=True requires xarray to be installed. Install with `pip install xarray`."
            ) from e
        return xr.Dataset(
            data_vars={
                "energy_rasters": (["iteration", "y", "x"], array_out[np.newaxis, 0]),
                "area_rasters": (["iteration", "y", "x"], array_out[np.newaxis, 1]),
                "energies": (["iteration"], np.array([self.energy], dtype=np.float64)),
                "temperatures": (["iteration"], np.array([temperature], dtype=np.float64)),
            },
            coords={
                "y": ("y", np.arange(self.dims[0])*self.resolution),
                "x": ("x", np.arange(self.dims[1])*self.resolution),
                "iteration": ("iteration", np.array([0], dtype=np.int32)),
            },
            attrs={
                "description": "OCN fit result arrays",
                "resolution_m": self.resolution,
                "gamma": self.gamma,
                "master_seed": int(self.master_seed),
            }
        )



    def fit(
        self,
        n_iterations:int=None,
        cooling_rate:float=1.0,
        constant_phase:float=0.0,
        pbar:bool=False,
        energy_reports:int=1000,
        array_reports:int=1,
        tol:float=None,
        max_iterations_per_loop=10_000,
        xarray_out:bool=False,
    ) -> "FitResult | xr.Dataset":
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
        energy_reports : int, default 1000
            Number of timepoints (approximately) at which energy is sampled. Must be > 0.
        array_reports : int, default 0
            Number of timepoints (approximately) at which the full grid is sampled to a numpy
            array. If None, do not report arrays. The returned array will have
            shape (array_reports, channels, rows, cols), where channel 0 (axis 1) is
            energy and channel 1 is drained area.
        tol : float, optional
            If provided, optimization will stop early if the relative change
            in energy between reports is less than `tol`. Must be positive.
            If None (default), no early stopping is performed.
            Recommended values are in the range 1e-4 to 1e-6.
        max_iterations_per_loop: int, optional
            If provided, the number of iterations steps to perform in each "chunk"
            of optimization. Energy and output arrays can be reported no more often
            than this. Recommended values are 1_000-1_000_000. Default is 10_000.
        xarray_out : bool, default False
            If True, the returned result will be an xarray.Dataset instead of a FitResult.
            Requires xarray to be installed.
        random_state : int | numpy.random.Generator | None, optional
            Seed or generator for RNG seeding. If None, system entropy is used.
            If an integer or Generator is provided, the new seed is drawn from it directly.

        Returns
        -------
        FitResult | tuple[FitResult, "xr.Dataset"]

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

        # validate inputs
        if n_iterations is None:
            n_iterations = int(40*self.dims[0]*self.dims[1])
        if not (isinstance(n_iterations, int) and n_iterations > 0):
            raise ValueError(f"n_iterations must be a positive integer, got {n_iterations}")
        n_iterations = max(energy_reports*10, n_iterations)

        if (not isinstance(tol, Number)) or tol < 0:
            if tol is not None:
                raise ValueError(f"tol must be a positive number or None, got {tol}")

        memory_est = array_reports*self.dims[0]*self.dims[1]*2*8 
        if memory_est > 20e6:
            warnings.warn(f"Requesting {array_reports} array is estimated to use {memory_est/1e6:.2f}MB of memory. Consider reducing array_reports or increasing max_iterations_per_loop if memory is a concern.")

        if (
            (cooling_rate < 0.5 or constant_phase > 0.1) 
            and n_iterations <= 50*self.dims[0]*self.dims[1]
        ):
            warnings.warn("Using cooling_rate < 0.5 and constant_phase > 0.1 with may cause convergence issues. Consider increasing n_iterations.")
        
        # make sure energy is up to date, useful if the user modified any parameters manually
        self.__p_c_graph.contents.energy = self.compute_energy()
        cooling_schedule = create_cooling_schedule(
            ocn=self,
            constant_phase=constant_phase,
            n_iterations=n_iterations,
            cooling_rate=cooling_rate,
        )

        # preallocate output arrays
        max_iterations_per_loop = int(max_iterations_per_loop)
        max_iterations_per_loop = max(1, max_iterations_per_loop)
        n_iterations = int(n_iterations)
        n_iterations = max(1, n_iterations)
        max_iterations_per_loop = min(n_iterations, max_iterations_per_loop)
        energy_report_interval = n_iterations // energy_reports
        energy_report_interval = max(energy_report_interval, max_iterations_per_loop)
        
        array_report_interval = n_iterations // array_reports
        array_report_interval = max(array_report_interval, max_iterations_per_loop)
        
        energy_out = np.empty(
            n_iterations//energy_report_interval + 2 + n_iterations//array_report_interval + 2, # always report energy when reporting arrays
            dtype=np.float64)
        energy_out[0] = self.energy

        array_out = np.empty([
            n_iterations//array_report_interval + 2, # adjust array_reports to be a multiple of n_loops
            2,
            self.dims[0],
            self.dims[1]
        ])
        array_out[0] = self.to_numpy()

        anneal_buf = np.empty(max_iterations_per_loop, dtype=np.float64)
        anneal_ptr = anneal_buf.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

        with tqdm(total=n_iterations, desc="OCN Optimization", unit_scale=True, dynamic_ncols=True, disable=not (pbar or self.verbosity >= 1)) as pbar:
            pbar.set_postfix({"Energy": self.energy, "P(Accept)": 1.0})
            pbar.update(0)
            
            array_report_number = 0
            energy_report_number = 0
            array_report_idx = [0]
            energy_report_idx = [0]
            
            self._advance_seed()
            
            completed_iterations = 0
            while completed_iterations < n_iterations:
                iterations_this_loop = min(max_iterations_per_loop, n_iterations - completed_iterations)
                anneal_buf[:iterations_this_loop] = cooling_schedule(
                    np.arange(completed_iterations, completed_iterations + iterations_this_loop)
                )

                e_old = self.energy
                check_status(_bindings.libocn.ocn_outer_ocn_loop(
                    self.__p_c_graph, 
                    iterations_this_loop, 
                    self.gamma, 
                    anneal_ptr
                ))
                e_new = self.energy
                completed_iterations += iterations_this_loop
                
                if (completed_iterations % array_report_interval) < max_iterations_per_loop:
                    array_report_number += 1
                    array_report_idx.append(completed_iterations)
                    array_out[array_report_number] = self.to_numpy()
                    energy_report_number += 1
                    energy_report_idx.append(completed_iterations)
                    energy_out[energy_report_number] = e_new
                elif (completed_iterations % energy_report_interval) < max_iterations_per_loop:
                    energy_report_number += 1
                    energy_report_idx.append(completed_iterations)
                    energy_out[energy_report_number] = e_new

                if tol is not None and e_new < e_old and abs((e_old - e_new)/e_old) < tol:
                    pbar.set_postfix({
                        "Energy": self.energy, 
                        "T": anneal_buf[iterations_this_loop - 1],
                        "Relative ΔE": (e_new - e_old)/e_old,
                    })
                    pbar.update(iterations_this_loop)
                    break 

                pbar.set_postfix({
                    "Energy": self.energy, 
                    "T": anneal_buf[iterations_this_loop - 1],
                    "Relative ΔE": (e_new - e_old)/e_old,
                })
                pbar.update(iterations_this_loop)

        array_report_idx = np.array(array_report_idx, dtype=np.int32)
        energy_report_idx = np.array(energy_report_idx, dtype=np.int32)
        
        if not xarray_out:
            return FitResult(
                i_raster=array_report_idx[:array_report_number + 1],
                energy_rasters=array_out[:array_report_number + 1, 0],
                area_rasters=array_out[:array_report_number + 1, 1],
                i_energy=energy_report_idx[:energy_report_number + 1],
                energies=energy_out[:energy_report_number + 1],
                temperatures=cooling_schedule(energy_report_idx[:energy_report_number + 1])
            )
        
        try:
            import xarray as xr
        except ImportError as e:
            raise ImportError(
                "PyOCN.OCN.fit() with xarray_out=True requires xarray to be installed. Install with `pip install xarray`."
            ) from e

        mask = np.isin(
            energy_report_idx[:energy_report_number + 1], 
            array_report_idx[:array_report_number + 1], 
            assume_unique=True
        )
        energy_idx = np.where(mask)[0]
        return xr.Dataset(
            data_vars={
                "energy_rasters": (
                    ["iteration", "y", "x"], 
                    array_out[:array_report_number + 1, 0]
                ),
                "area_rasters": (
                    ["iteration", "y", "x"], 
                    array_out[:array_report_number + 1, 1]
                ),
                "energies": (["iteration"], energy_out[energy_idx]),
                "temperatures": (["iteration"], cooling_schedule(array_report_idx[:array_report_number + 1]))
            },
            coords={
                "iteration": ("iteration", array_report_idx),
                "y": ("y", np.arange(self.dims[0])*self.resolution),
                "x": ("x", np.arange(self.dims[1])*self.resolution),
            },
            attrs={
                "description": "OCN fit result arrays",
                "resolution_m": self.resolution,
                "gamma": self.gamma,
                "master_seed": int(self.master_seed),
            }
        )


__all__ = [
    "OCN",
]
