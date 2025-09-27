from __future__ import annotations
from itertools import product
from typing import Literal, Callable, TYPE_CHECKING
import networkx as nx
import numpy as np
from tqdm import tqdm

if TYPE_CHECKING:
    from .ocn import OCN

_allowed_net_types = {"I", "H", "V", "T"}

#TODO: add ability to move root?
def net_type_to_dag(net_type:Literal["I", "H", "V", "T"], dims:tuple, pbar=False) -> nx.DiGraph:
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
        case "T":
            raise NotImplementedError("T net type not yet implemented.")
        case _:
            raise ValueError(f"Invalid net_type {net_type}. Must be one of {_allowed_net_types}.")
        
    return G

def create_cooling_schedule(
    ocn:OCN,
    constant_phase:float, 
    n_iterations:int, 
    cooling_rate:float, 
) -> Callable[[int], float|np.ndarray]:
    """
    Simulated annealing temperature. 
    The function that expresses the temperature of the simulated annealing process is as follows:

    if i <= initialNoCoolingPhase*nIter:
    Temperature[i] = Energy[1]

    if initialNoCoolingPhase*nIter < i <= nIter:
    Temperature[i] = Energy[1]*(-coolingRate*(i - InitialNocoolingPhase*nIter)/nNodes)

    where i is the index of the current iteration 
    and Energy[1] = sum(A^expEnergy), 
    with A denoting the vector of drainage areas corresponding to the initial state of the network. 
    According to the simulated annealing principle, a new network configuration obtained at iteration i 
    is accepted with probability equal to exp((Energy[i] - Energy[i-1])/Temperature[i]) 
    if Energy[i] < Energy[i-1]. 
    To ensure convergence, it is recommended to use coolingRate values between 0.5 and 10 and initialNoCoolingPhase <= 0.3. 
    Low coolingRate and high initialNoCoolingPhase values cause the network configuration to depart more significantly from the initial state. 
    If coolingRate < 0.5 and initialNoCoolingPhase > 0.1 are used, it is suggested to increase nIter with respect to the default value in order to guarantee convergence.
    """
    n_constant = int(constant_phase * n_iterations)
    nnodes = ocn.dims[0] * ocn.dims[1]

    term1 = -cooling_rate / nnodes
    term2 = cooling_rate * n_constant / nnodes

    e0 = ocn.energy

    def schedule(i):
        return np.where(i < n_constant, e0, e0 * np.exp(term1*i + term2))

    return schedule