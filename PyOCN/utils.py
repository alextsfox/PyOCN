"""
Utility functions for working with OCNs.
"""

from __future__ import annotations
from itertools import product
from typing import Literal, Callable, TYPE_CHECKING
import networkx as nx
import numpy as np
from tqdm import tqdm

if TYPE_CHECKING:
    from .ocn import OCN

_allowed_net_types = {"I", "H", "V", "T", "E"}

#TODO: add ability to move root?
def net_type_to_dag(net_type:Literal["I", "H", "V", "T", "E"], dims:tuple, pbar: bool = False) -> nx.DiGraph:
    """Create a predefined OCN initialization network as a NetworkX DiGraph.

    Parameters
    ----------
    net_type : {"I", "H", "V", "T", "E"}
        The type of network to create.
        Descriptions of allowed types:

        - "I":

          ::

              O--O--O--O--O
                    |
              O--O--O--O--O
                    |
              O--O--O--O--O
                    |
              O--O--X--O--O

        - "V":

          ::

              O  O  O  O  O
               \  \ | /  /
              O  O  O  O  O
               \  \ | /  /
              O  O  O  O  O
               \  \ | /  /
              O--O--X--O--O

        - "H":

          ::

              O  O  O  O
              |  |  | /
              O  O  O--O
              |  | /
              O  O--O--O
              | /
              X--O--O--O

        - "E": A network where every node on the edge of the grid is a root.

          ::
              X  X  X  X  X  X
                \  \    /  /
              X  O  O  O  O  X
                \  \    /  /
              X  O  O  O  O  X
  
              X  O  O  O  O  X
                /  /     \  \
              X  O  O  O  O  X
                /  /     \  \  
              X  X  X  X  X  X


        - "T": Not implemented yet.

    dims : tuple
        The network dimensions as ``(rows, cols)``. Both must be positive even integers.
    pbar : bool, default False
        If True, display a progress bar while constructing the graph.

    Returns
    -------
    networkx.DiGraph
        A directed acyclic graph representing a valid initial OCN configuration.

    Raises
    ------
    ValueError
        If ``net_type`` is invalid or ``dims`` are not two positive even integers.

    Notes
    -----
    The returned graph assigns each grid cell exactly one node with a ``pos``
    attribute equal to ``(row, col)``.
    """
    rows, cols = dims
    G = nx.DiGraph()
    pbar = tqdm(product(range(rows), range(cols)), total=rows*cols, disable=not pbar, desc="Creating DAG")
    match net_type:
        case "I":
            jroot = cols // 2
            for i, j in pbar:
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
            for i, j in pbar:
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
            for i, j in pbar:
                n = i*cols + j
                G.add_node(n, pos=(i, j))
                if i == j and i > 0:  # main diagonal
                    G.add_edge(n, n - cols - 1)
                elif i > j:
                    G.add_edge(n, n - cols)
                elif j > i:
                    G.add_edge(n, n - 1)
        case "E":  #TODO: implement a better radial pattern
            half_rows = rows / 2
            half_cols = cols / 2
            for i, j in pbar:
                n = i*cols + j
                G.add_node(n, pos=(i, j))
                # ul quadrant: flow up and left
                if i < half_rows and j < half_cols and i > 0 and j > 0:
                    G.add_edge(n, n - cols - 1)
                # ur quadrant: flow up and right
                elif i < half_rows and j >= half_cols and i > 0 and j < cols - 1:
                    G.add_edge(n, n - cols + 1)
                # ll quadrant: flow down and left
                elif i >= half_rows and j < half_cols and j > 0 and i < rows - 1:
                    G.add_edge(n, n + cols - 1)
                # lr quadrant: flow down and right
                elif i >= half_rows and j >= half_cols and i < rows - 1 and j < cols - 1:
                    G.add_edge(n, n + cols + 1)
        case _:
            raise ValueError(f"Invalid net_type {net_type}. Must be one of {_allowed_net_types}.")
        
    return G

def create_cooling_schedule(
    ocn: OCN,
    constant_phase: float,
    n_iterations: int,
    cooling_rate: float,
) -> Callable[[int], float | np.ndarray]:
    """Create a simulated-annealing cooling schedule for OCN optimization.

    This returns a callable ``schedule(i)`` that returns the temperature at
    iteration ``i``. The schedule consists of a constant-temperature phase
    followed by an exponentially decaying phase.

    Parameters
    ----------
    ocn : OCN
        The OCN instance (used to infer the number of nodes, ``rows*cols``).
    constant_phase : float
        Fraction of iterations (``0 <= fraction <= 1``) during which the
        temperature remains constant at ``Energy[0]``.
    n_iterations : int
        Total number of optimization iterations.
    cooling_rate : float
        Positive decay rate controlling the exponential temperature decrease
        after the constant phase.

    Returns
    -------
    Callable[[int], float] | numpy.ndarray
        A function mapping an iteration index ``i`` to a temperature value. If
        vectorized evaluation is used, may return a NumPy array of temperatures.

    Notes
    -----
    The exponential phase follows the form

    .. math::

        T_i = E_0 \exp\left(-\frac{\text{cooling\_rate}\,(i - n_0)}{N}\right),

    where ``E0`` is the initial energy, ``n0`` is the number of iterations in
    the constant phase, and ``N = rows * cols``.
    """
    n_constant = int(constant_phase * n_iterations)
    nnodes = ocn.dims[0] * ocn.dims[1]

    term1 = -cooling_rate / nnodes
    term2 = cooling_rate * n_constant / nnodes

    e0 = ocn.energy

    def schedule(i):
        return np.where(i < n_constant, e0, e0 * np.exp(term1*i + term2))

    return schedule