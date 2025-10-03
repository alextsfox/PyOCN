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

"""
TODO de-obfuscate the return values of the fit methods. 
I think include the external fit functions as methods of OCN
The fit methods always save the energy to an attribute,
and optionally return arrays if reports are requested.
"""

import warnings
import ctypes
from typing import Any, Callable
from os import PathLike
from numbers import Number
from pathlib import Path
from dataclasses import dataclass

import networkx as nx 
import numpy as np
from tqdm import tqdm

from ._statushandler import check_status
from .utils import simulated_annealing_schedule, net_type_to_dag, unwrap_digraph, assign_subwatersheds
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
    watershed_id : np.ndarray
        Array snapshots of watershed IDs assigned during fitting, shaped
        (n_reports, rows, cols). Roots are assigned -1. Non-draining cells are -9999.
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
    watershed_id: np.ndarray
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
        Export raster arrays (energy, drained area, watershed_id) as numpy arrays.
    :meth:`to_xarray`
        Export raster arrays as an xarray Dataset (requires xarray).
    :meth:`to_gtiff`
        Export raster arrays to a GeoTIFF file (requires rasterio).
    :meth:`copy`
        Create a deep copy of the OCN.

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
    wrap : bool
        If true, allows wrapping around the edges of the grid (toroidal). If false, no wrapping is applied (read-only property).
    history : np.ndarray
        numpy array of shape (n_iterations, 3) recording the iteration index, energy, and temperature at each iteration during optimization.
        If multiple fit calls are made, history is appended to this array.

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
    >>> TODO: update this example to match new fit API
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
    def __init__(self, dag: nx.DiGraph, resolution: float=1.0, gamma: float = 0.5, random_state=None, verbosity: int = 0, validate:bool=True, wrap : bool = False):
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
        # sets nroots, resolution, dims, wrap
        self.__p_c_graph = fgconv.from_digraph(dag, resolution, verbose=(verbosity > 1), validate=validate, wrap=wrap)
        self.__p_c_graph.contents.energy = self.compute_energy()

        self.__history = np.empty((0, 3), dtype=np.float64)

    @classmethod
    def from_net_type(cls, net_type:str, dims:tuple, resolution:float=1, gamma : float = 0.5, random_state=None, verbosity:int=0, wrap : bool = False):
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
        wrap : bool, default False
            If true, allows wrapping around the edges of the grid (toroidal). If false, no wrapping is applied.

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
        return cls(dag, resolution, gamma, random_state, verbosity=verbosity, validate=False, wrap=wrap)

    @classmethod
    def from_digraph(cls, dag: nx.DiGraph, resolution:float=1, gamma=0.5, random_state=None, verbosity: int = 0, wrap: bool = False):
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
        wrap : bool, default False
            If true, allows wrapping around the edges of the grid (toroidal). If false, no wrapping is applied.

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
          i.e., no jumps over rows or columns. If ``wrap=True``, edges may
          connect across the grid boundaries (i.e. row 0 can connect to row m-1 and col 0 can connect to col n-1).
        - Edges do not cross in the row-column plane.
        - Both ``m`` and ``n`` are even integers, and there are at least four
          vertices.
        """
        return cls(dag, resolution, gamma, random_state, verbosity=verbosity, validate=True, wrap=wrap)

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
            self.__history.nbytes +
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
        cpy.__history = self.history.copy()
        cpy.__p_c_graph.contents.resolution = self.resolution
        cpy.__p_c_graph.contents.wrap = self.wrap

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
    
    def copy(self) -> "OCN":
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
    def wrap(self) -> bool:
        """
        Whether the grid allows wrapping around the edges (toroidal).

        Returns
        -------
        bool
            Current wrap setting.
        """
        return self.__p_c_graph.contents.wrap
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
        if not isinstance(random_state, (int, np.integer, type(None), np.random.Generator)):
            raise ValueError("RNG must be initialized with an integer/Generator/None.")
        self.__master_seed = np.random.default_rng(random_state).integers(0, int(2**32 - 1))
        _bindings.libocn.rng_seed(self.master_seed)
    def _advance_seed(self):
        random_state = self.master_seed
        new_random_state = np.random.default_rng(random_state).integers(0, int(2**32 - 1))
        self.__master_seed = new_random_state
        _bindings.libocn.rng_seed(self.__master_seed)
    @property
    def history(self) -> np.ndarray:
        """
        numpy array of shape (n_iterations, 3) recording the iteration index, energy, and temperature at each iteration during optimization.
        If multiple fit calls are made, history is appended to this array.

        Returns
        -------
        np.ndarray
            The optimization history.
        """
        return self.__history

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
            - ``watershed_id``: integer watershed ID (roots have watershed id = -1)
        """
        dag = fgconv.to_digraph(self.__p_c_graph.contents)
        assign_subwatersheds(dag)

        node_energies = dict()

        for node in nx.topological_sort(dag):
            node_energies[node] = (
                dag.nodes[node]['drained_area']**self.gamma 
                + sum(node_energies[p] for p in dag.predecessors(node))
            )
        nx.set_node_attributes(dag, node_energies, 'energy')

        return dag
    
    def to_gtiff(self, west:float, north:float, crs: Any, path:str|PathLike, unwrap:bool=True):
        """
        Export a raster of the current FlowGrid to a .gtiff file
        using rasterio. The resulting raster has 3 bands: `energy`, `drained_area`, and `watershed_id`.
        The `watershed_id` band contains integer watershed IDs assigned to each node,
        with root nodes assigned a value of -1. NA values are either np.nan (for energy and drained_area)
        or -9999 (for watershed_id).

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
            The northern border of the raster in crs units, corresponding
            to row 0.
        crs : Any
            The crs for the resulting gtiff, passed to `rasterio.open`
        path : str or Pathlike
            The output path for the resulting gtiff file.
        unwrap : bool, default True
            If True, unwraps the DAG to a continuous dag with non-periodic
            boundaries before exporting. This will result in a larger raster if
            the current OCN is wrapping.
        """
        try:
            import rasterio
            from rasterio.transform import from_origin
        except ImportError as e:
            raise ImportError(
                "PyOCN.OCN.to_gtiff() requires rasterio to be installed. Install with `pip install rasterio`."
            ) from e

        array = self.to_numpy(unwrap=unwrap)
        dims = array.shape[1:]
        energy = array[0]
        drained_area = array[1]
        watershed_id = np.where(np.isnan(array[2]), -9999, array[2])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            watershed_id = watershed_id.astype(np.int32)

        transform = from_origin(west, north, self.resolution, self.resolution)

        # Write three bands: 1=energy, 2=drained_area, 3=watershed_id
        with rasterio.open(
            Path(path),
            "w",
            driver="GTiff",
            height=dims[0],
            width=dims[1],
            count=3,
            dtype=np.float64,
            crs=crs,
            transform=transform,
            compress="deflate",
        ) as dst:
            dst.write(energy, 1)
            dst.write(drained_area, 2)
            dst.write(watershed_id, 3)
            # Band descriptions (nice to have)
            try:
                dst.set_band_description(1, "energy")
                dst.set_band_description(2, "drained_area")
                dst.set_band_description(3, "watershed_id")
            except Exception:
                pass
    
    def to_numpy(self, unwrap:bool=True) -> np.ndarray:
        """
        Export the current FlowGrid to a numpy array with shape (2, rows, cols).
        Has two channels: 0=energy, 1=drained_area.

        Parameters
        ----------
        unwrap : bool, default True
            If True, unwraps the DAG to a non-wrapping representation
            before exporting. This will result in a larger raster if
            the current OCN is wrapping.
        """
        dag = self.to_digraph()
        dims = self.dims
        if self.wrap and unwrap:
            dag = unwrap_digraph(dag, dims)
            positions = np.array(list(nx.get_node_attributes(dag, 'pos').values()))
            max_r, max_c = positions.max(axis=0)
            dims = (max_r + 1, max_c + 1)

        energy = np.full(dims, np.nan)
        drained_area = np.full(dims, np.nan)
        watershed_id = np.full(dims, np.nan)
        for node in dag.nodes:
            r, c = dag.nodes[node]['pos']
            energy[r, c] = dag.nodes[node]['energy']
            drained_area[r, c] = dag.nodes[node]['drained_area']
            watershed_id[r, c] = dag.nodes[node]['watershed_id']
        return np.stack([energy, drained_area, watershed_id], axis=0)

    def to_xarray(self, unwrap:bool=True) -> "xr.Dataset":
        """
        Export the current FlowGrid to an xarray Dataset
        
        Parameters
        ----------
        unwrap : bool, default True
            If True, unwraps the DAG to a non-wrapping representation
            before exporting. This will result in a larger raster if
            the current OCN is wrapping. If unwrapping, the (0,0) coordinate
            will be set to the position of the "main" root node, defined as
            the root node with the smallest row*cols + col value. Otherwise,
            (0,0) will be the top-left corner of the grid.
        
        Returns
        -------
        xr.Dataset
         an xarray Dataset with data variables:
            - `energy_rasters` (np.float64) representing energy at each grid cell
            - `area_rasters` (np.float64) representing drained area at each grid cell
            - `watershed_id` (np.int32). NA value is -9999. Roots have value -1. Represents the watershed membership ID for each grid cell.
        and coordinates:
            - `y` (float) representing the northing coordinate of each row in meters
            - `x` (float) representing the easting coordinate of each column in meters
        """

        try:
            import xarray as xr
        except ImportError as e:
            raise ImportError(
                "PyOCN.OCN.to_xarray() requires xarray to be installed. Install with `pip install xarray`."
            ) from e
        
        dims = self.dims

        dag = self.to_digraph()
        row_root, col_root = 0, 0
        if self.wrap and unwrap:
            roots = [n for n, d in dag.out_degree() if d==0]
            main_root = min(roots, key=lambda n: dag.nodes[n]['pos'][0]*dims[1] + dag.nodes[n]['pos'][1])
            
            dag = unwrap_digraph(dag, dims)
            positions = np.array(list(nx.get_node_attributes(dag, 'pos').values()))
            max_r, max_c = positions.max(axis=0)
            dims = (max_r + 1, max_c + 1)
            
            # compute the new position of the root node after unwrapping. This will be the new origin (0,0).
            row_root, col_root = dag.nodes[main_root]['pos']

        energy = np.full(dims, np.nan)
        drained_area = np.full(dims, np.nan)
        watershed_id = np.full(dims, np.nan)
        for node in dag.nodes:
            r, c = dag.nodes[node]['pos']
            energy[r, c] = dag.nodes[node]['energy']
            drained_area[r, c] = dag.nodes[node]['drained_area']
            watershed_id[r, c] = dag.nodes[node]['watershed_id']
        array_out = np.stack([energy, drained_area, watershed_id], axis=0)

        dims = array_out.shape[1:]

        watershed_id = np.where(np.isnan(array_out[2]), -9999, array_out[2])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            watershed_id = watershed_id.astype(np.int32)

        return xr.Dataset(
            data_vars={
                "energy_rasters": (["y", "x"], array_out[0].astype(np.float64)),
                "area_rasters": (["y", "x"], array_out[1].astype(np.float64)),
                "watershed_id": (["y", "x"], watershed_id),
            },
            coords={
                "y": ("y", np.linspace(-row_root, (-row_root + (dims[0]-1)), dims[0])*self.resolution),
                "x": ("x", np.linspace(-col_root, (-col_root + (dims[1]-1)), dims[1])*self.resolution),
            },
            attrs={
                "description": "OCN fit result arrays",
                "resolution_m": self.resolution,
                "gamma": self.gamma,
                "master_seed": int(self.master_seed),
                "wrap": self.wrap,
            }
        )

    def single_iteration(self, temperature:float, array_report:bool=False) -> "xr.Dataset | None":
        """ 
        Perform a single iteration of the optimization algorithm at a given temperature. Updates the internal history attribute.
        See :meth:`fit` for details on the algorithm.

        Parameters
        ----------
        temperature : float
            Temperature parameter governing acceptance probability. Typical
            range is a fraction of ocn.energy.
        array_report : bool, default False
            If True (default), the returned result will be an xarray.Dataset.
            See :meth:`fit` for details. The returned object will have an iteration dimension of size 1.
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
            temperature,
        ))

        history = np.array([[self.__history[-1, 0] + 1, self.energy, temperature]])
        self.__history = np.concatenate((self.__history, history), axis=0)
        
        if not array_report:
            return None
        
        ds = self.to_xarray(unwrap=self.unwrap)
        ds = ds.expand_dims({"iteration": [0]})

        return ds
    
    def fit(
        self,
        cooling_rate:float=1.0,
        constant_phase:float=0.0,
        n_iterations:int=None,
        pbar:bool=False,
        array_reports:int=0,
        tol:float=None,
        max_iterations_per_loop=10_000,
    ) -> "xr.Dataset | None":
        """
        Convenience function to optimize the OCN using the simulated annealing algorithm from Carraro et al (2020).
        For finer control over the optimization process, use :meth:`fit_custom_cooling` or use :meth:`single_erosion_event` in a loop.

        This performs ``n_iterations`` erosion events, accepting or rejecting
        proposals according to a temperature schedule defined by the annealing algorithm. 
        A proposal consists of changing the outflow direction of a randomly selected vertex. 
        The new outflow direction is chosen uniformly from the valid neighbors.
        A proposal is valid if it maintains a well-formed graph structure.

        The :attr:`history` attribute is updated in-place after optimization finishes.

        Parameters
        ----------
        ocn : OCN
            The OCN instance to optimize.
        cooling_rate : float, default 1.0
            Cooling rate parameter in the annealing algorithm. Typical range is 0.5-1.5.
            Higher values result in faster cooling and a greedier search.
            Lower values result in slower cooling and more thorough exploration of the solution space, but slower convergence and lower stability.
        constant_phase : float, default 0.0
            Amount of time to hold temeprature constant at the start of the optimization.
            This is a fraction of n_iterations, and must be in the range [0, 1].
            A value of 0.0 (default) means the temperature starts cooling immediately
            from the initial temperature. A value of 1.0 means the temperature is held
            constant for the entire optimization.
        n_iterations : int, optional
            Total number of iterations. Defaults to ``40 * rows * cols``.
            Always at least ``energy_reports * 10`` (this should only matter for
            extremely small grids, where ``rows * cols < 256``).
        pbar : bool, default True
            Whether to display a progress bar.
        array_reports : int, default 0
            Number of timepoints (approximately) at which to save the state of the FlowGrid.
            If 0 (default), returns None. If >0, returns a FitResult or xarray.Dataset
            containing the state of the FlowGrid at approximately evenly spaced intervals
            throughout the optimization, including the initial and final states.
            
            array_reports > 0, the returned result will be an xarray.Dataset. Requires xarray to be installed.
            
            The returned xarray.Dataset will have coordinates:
                - `y` (float) representing the northing coordinate of each row in meters
                - `x` (float) representing the easting coordinate of each column in meters
                - `iteration` (int) representing the iteration index at which the data was recorded
            data variables:
                - `energy_rasters` (np.float64) representing energy at each grid cell
                - `area_rasters` (np.float64) representing drained area at each grid cell
                - `watershed_id` (np.int32). NA value is -9999. Roots have value -1. Represents the watershed membership ID for each grid cell.
            The coordiante (0, 0) is the top-left corner of the grid.

            If the OCN has a periodic boundary condition, the following changes apply: 
                - The (0,0) coordinate will be set to the position of the "main" root node, defined as
                the root node with the smallest row*cols + col value
                - The rasters will be unwrapped to a non-periodic representation, which may result in larger rasters.
                - The size of the final rasters are the maximum extent of the unwrapped grid, taken across all iterations.

            Generating reports requires additional memory and computation time.
            
        tol : float, optional
            If provided, optimization will stop early if the relative change
            in energy between reports is less than `tol`. Must be positive.
            If None (default), no early stopping is performed.
            Recommended values are in the range 1e-4 to 1e-6.
        max_iterations_per_loop: int, optional
            If provided, the number of iterations steps to perform in each "chunk"
            of optimization. Energy and output arrays can be reported no more often
            than this. Recommended values are 1_000-1_000_000. Default is 10_000.

        Returns
        -------
        FitResult | xr.Dataset | None

        Raises
        ------
        ValueError

        Warns
        -----
        UserWarning

        Optimization Algorithm
        ----------------------
        At iteration ``i``, the outflow of a random grid cell if proposed to be rerouted.
        The proposal is accepted with the probability
        .. math::

            P(\text{accept}) = e^{-\Delta E / T},

        where :math:`\Delta E` is the change in energy the change would cause 
        and :math:`T` is the temperature of the network.

        The total energy of the system is computed from the drained areas of each grid cell :math:`k` as

        .. math::
            E = \sum_k A_k^\gamma

        The temperature of the network is governed by a cooling schedule, which is a function of iteration index.
        
        Note that when :math:`\Delta E < 0`, the move is always accepted.

        Simulated Annealing Schedule
        ----------------------------
        The cooling schedule used by this method is a piecewise function of iteration index:
        .. math::

            T(i) = \begin{cases}
                E_0 & i < C N \\
                E_0 \cdot e^{\;i - C N} & i \ge C N
            \end{cases}

        where :math:`E_0` is the initial energy, :math:`N` is the total number
        of iterations, and :math:`C` is ``constant_phase``. Decreasing-energy
        moves (:math:`\Delta E < 0`) are always accepted.

        Alternative cooling schedules can be implemented using :meth:`fit_custom_cooling`.
        """
        
        # make sure energy is up to date, useful if the user modified any parameters manually
        self.__p_c_graph.contents.energy = self.compute_energy()
        if constant_phase is None:
            constant_phase = 0.0
        if cooling_rate is None:
            cooling_rate = 1.0
        if n_iterations is None:
            n_iterations = 40 * self.dims[0] * self.dims[1]
        cooling_func = simulated_annealing_schedule(
            dims=self.dims,
            E0=self.energy,
            constant_phase=constant_phase,
            n_iterations=n_iterations,
            cooling_rate=cooling_rate,
        )

        return self.fit_custom_cooling(
            cooling_func=cooling_func,
            n_iterations=n_iterations,
            pbar=pbar,
            array_reports=array_reports,
            tol=tol,
            max_iterations_per_loop=max_iterations_per_loop,
        )

    def fit_custom_cooling(
        self,
        cooling_func:Callable[[np.ndarray], np.ndarray],
        n_iterations:int=None,
        iteration_start:int=0,
        pbar:bool=False,
        array_reports:int=0,
        tol:float=None,
        max_iterations_per_loop=10_000,
    ) -> "xr.Dataset | None":
        """
        Optimize the OCN using the a custom cooling schedule. This allows for
        multi-stage optimizations or other custom cooling schedules not covered by the default simulated annealing schedule
        from Carraro et al (2020).

        See :meth:`fit` for additional details on the optimization algorithm and parameters.

        Parameters
        ----------
        ocn : OCN
            The OCN instance to optimize.
        cooling_func : Callable[[np.ndarray], np.ndarray]
            A function that takes an array of iteration indices and returns an array of temperatures.
            This function defines the cooling schedule for the optimization. Note that the function
            should return temperatures that are appropriate for the current energy of the OCN.
        n_iterations : int, optional
        iteration_start : int, default 0
            The starting iteration index. This is useful for continuing an optimization
            from a previous run. Must be a non-negative integer. If provided, ``n_iterations``
            is the number of additional iterations to perform. The iteration number passed to
            ``cooling_func`` will be ``iteration_start + i`` where ``i`` is the current iteration index
            in the range ``[0, n_iterations-1]``.
        pbar : bool, default True
        array_reports : int, default 0
        tol : float, optional
        max_iterations_per_loop: int, optional

        Returns
        -------
        FitResult | xr.Dataset | None

        Raises
        ------
        ValueError

        Warns
        -----
        UserWarning
        """
        # validate inputs
        if n_iterations is None:
            n_iterations = int(40*self.dims[0]*self.dims[1])
        if not (isinstance(n_iterations, int) and n_iterations > 0):
            raise ValueError(f"n_iterations must be a positive integer, got {n_iterations}")
        if not (isinstance(array_reports, int) and array_reports >= 0):
            raise ValueError(f"array_reports must be a non-negative integer, got {array_reports}")
        if (not isinstance(iteration_start, (int, np.integer))) or iteration_start < 0:
            raise ValueError(f"iteration_start must be a non-negative integer, got {iteration_start}")
        if (not isinstance(tol, Number)) or tol < 0:
            if tol is not None:
                raise ValueError(f"tol must be a positive number or None, got {tol}")

        xarray_out = array_reports > 0

        memory_est = array_reports*self.dims[0]*self.dims[1]*2*8 
        if memory_est > 20e6:
            warnings.warn(f"Requesting {array_reports} array is estimated to use {memory_est/1e6:.2f}MB of memory. Consider reducing array_reports or increasing max_iterations_per_loop if memory is a concern.")
        
        self.__p_c_graph.contents.energy = self.compute_energy()
        cooling_schedule = cooling_func

        # preallocate output arrays
        max_iterations_per_loop = int(max_iterations_per_loop)
        max_iterations_per_loop = max(1, max_iterations_per_loop)
        n_iterations = int(n_iterations)
        n_iterations = max(1, n_iterations)
        max_iterations_per_loop = min(n_iterations, max_iterations_per_loop)
        

        # set up output arrays
        energy_out = np.empty(
            n_iterations//max_iterations_per_loop + 2, # always report energy when reporting arrays
            dtype=np.float64)
        
        if array_reports:
            array_report_interval = n_iterations // array_reports
            array_report_interval = max(array_report_interval, max_iterations_per_loop)
        else:
            array_report_interval = n_iterations*2  # never report arrays
        if array_reports and not xarray_out:  # we use a different method for xarray output
            array_out = np.empty([
                n_iterations//array_report_interval + 2, # adjust array_reports to be a multiple of n_loops
                3,
                self.dims[0],
                self.dims[1]
            ])
            array_out[0] = self.to_numpy(unwrap=False)

        anneal_buf = np.empty(max_iterations_per_loop, dtype=np.float64)
        anneal_ptr = anneal_buf.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

        with tqdm(total=n_iterations, desc="OCN Optimization", unit_scale=True, dynamic_ncols=True, disable=not (pbar or self.verbosity >= 1)) as pbar:
            pbar.set_postfix({"Energy": self.energy, "P(Accept)": 1.0})
            pbar.update(0)
            
            self._advance_seed()

            completed_iterations = iteration_start
            n_iterations += iteration_start

            # set up reporting
            array_report_idx = []
            energy_report_idx = []
            ds_out_dict = dict()
            array_report_idx.append(completed_iterations)
            ds_out_dict[completed_iterations] = self.to_xarray(unwrap=self.wrap)
            energy_out[0] = self.energy
            energy_report_idx.append(completed_iterations)
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
                    anneal_ptr,
                ))
                e_new = self.energy
                completed_iterations += iterations_this_loop
                
                # report arrays (and energy)
                if (
                    xarray_out
                    and (
                        (completed_iterations % array_report_interval) < max_iterations_per_loop
                        or completed_iterations >= n_iterations
                    )
                ):
                    array_report_idx.append(completed_iterations)
                    ds_out_dict[completed_iterations] = self.to_xarray(unwrap=self.wrap)

                # always report energy
                energy_out[len(energy_report_idx)] = e_new
                energy_report_idx.append(completed_iterations)

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

        # save energy and temp history
        start = 1 if self.__history.shape[0] > 0 else 0  # avoid duplicating the initial state if continuing from previous fit
        idx_offset = self.__history[-1, 0] if self.__history.shape[0] > 0 else 0
        history = np.stack([
            np.array(idx_offset - iteration_start + np.asarray(energy_report_idx)[start:]),
            np.array(energy_out[start:len(energy_report_idx)]),
            np.array(cooling_schedule(energy_report_idx[start:]))
        ], axis=1)
        self.__history = np.concatenate([self.__history, history], axis=0)
        
        if not xarray_out: 
            return
    
        try:
            import xarray as xr
        except ImportError as e:
            raise ImportError(
                "PyOCN.OCN.fit() with array_report>0 requires xarray to be installed. Install with `pip install xarray`."
            ) from e
        
        # find the maximum extent of the unwrapped grid across all reported arrays
        coord_ranges = list((ds.x.data.min(), ds.x.data.max(), ds.y.data.min(), ds.y.data.max()) for ds in ds_out_dict.values())
        xmin, xmax, ymin, ymax = (
            min(coord_range[0] for coord_range in coord_ranges), 
            max(coord_range[1] for coord_range in coord_ranges), 
            min(coord_range[2] for coord_range in coord_ranges), 
            max(coord_range[3] for coord_range in coord_ranges),
        )
        new_xcoords = np.arange(xmin, xmax + self.resolution, self.resolution)
        new_ycoords = np.arange(ymin, ymax + self.resolution, self.resolution)

        data_shape = (len(ds_out_dict), len(new_ycoords), len(new_xcoords))


        # build an empty dataset
        ds = xr.Dataset(
            data_vars={
                "energy_rasters": (
                    ["iteration", "y", "x"], 
                    np.full(data_shape, np.nan, dtype=np.float64)
                ),
                "area_rasters": (
                    ["iteration", "y", "x"], 
                    np.full(data_shape, np.nan, dtype=np.float64)
                ),
                "watershed_id": (
                    ["iteration", "y", "x"], 
                    np.full(data_shape, -9999, dtype=np.int32)
                ),
            },
            coords={
                "iteration": ("iteration", np.asarray(array_report_idx)),
                "y": ("y", new_ycoords),
                "x": ("x", new_xcoords),
            },
            attrs={
                "description": "OCN fit result arrays",
                "resolution_m": self.resolution,
                "gamma": self.gamma,
                "master_seed": int(self.master_seed),
                "wrap": self.wrap,
            }
        )

        # fill in the dataset with the unwrapped arrays
        for i, ds_i in ds_out_dict.items():
            ds.energy_rasters.loc[dict(iteration=i, y=ds_i.y, x=ds_i.x)] = ds_i.energy_rasters
            ds.area_rasters.loc[dict(iteration=i, y=ds_i.y, x=ds_i.x)] = ds_i.area_rasters
            ds.watershed_id.loc[dict(iteration=i, y=ds_i.y, x=ds_i.x)] = ds_i.watershed_id

        return ds
