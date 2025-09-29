# PyOCN
This is a package to generate optimal channel networks (OCNs), based on the algorithm described in Carraro et al. (2020). *Generation and application of river network analogues for use in ecology and evolution. Ecology and Evolution.* doi:10.1002/ece3.6479 and mirrors some of the functionality of the OCNet R package (https://lucarraro.github.io/OCNet/).

This is a work in progress. I plan to release this as a package on PyPI in the near future. For now, you can clone the repository and build the package from source.

To compile libocn with `gcc` on MacOS or Linux, move the root directory into your project's working directory andrun the following commands from the root directory:

```bash
cd PyOCN/c_src
bash build.sh
```

or alternatively:

```bash
cd PyOCN/c_src
gcc -fPIC -O3 -flto -c ocn.c flowgrid.c status.c rng.c
gcc -shared -O3 -flto -o libocn.so ocn.o flowgrid.o status.o rng.o
mv libocn.so ../libocn.so
rm -f ocn.o flowgrid.o status.o rng.o
```

If you have any questions or comments, please open an issue or contact me directly at https://www.afox.land

# OCN Algorithm
An initial stream network is generated as a directed acyclic graph (DAG) that is a spanning tree over a 2d grid of cells. Each cell in the grid (except the root) has a single outgoing edge that connects to one of its 8 neighboring cells. The OCN algorithm then iteratively modifies the DAG by randomly selecting a cell and changing its outflow to point to a different neighbor, ensuring that the spanning tree structure is maintained (*i.e. no cycles are introduced and each cell has a single outflow, except for the root cell*). The total energy at cell $k$ after iteration $n$ is given by:

$$
E_k[n] = \sum_{i} A_i^\gamma[n]
$$

Where $A_i$ is the cumulative drained area at cell $i$, and $i\in\mathrm{drained}(k)$ is a cell drained by $k$, including itself, and $\gamma$ controls how the energy scales with area. This change is accepted or rejected based on an annealing method, which is related to the Metropolis-Hastings algorithm. The probability of accepting a proposed change to the network is given by:

$$
P(\mathrm{accept}) = \min\left(1, \exp\left(-\frac{E_\mathrm{root}[n] - E_\mathrm{root}[n-1]}{T[n]}\right)\right)
$$

Where $T[n]$ is the "temperature" of the network at iteration $n$. Initially, the temperature is set to a high value ($\sim E[0]$) to encourage exploration of the solution space. The temperature is then set to exponentially decay over time, "annealing" the network as it settles into a low-energy configuration. Note that a proposed change is always accepted if it results in a lower energy state ($E[n] < E[n-1]$).

Lower values of gamma allow very dendritic networks with lots fo branching. The following animation shows the process of optimizing an OCN with $\gamma=0.7$ on a 256x256 grid for 4M iterations, which took about 2 minutes to run on a Macbook pro. 

<div align="center">
  <img src="generation.gif" alt="Optimizing an OCN">
</div>

# libocn
The backend of PyOCN is the libocn C library. libocn implements the core algorithms for generating and manipulating OCNs. Unlike the OCNet R package which uses an adjacency matrix implementation representation of the flow grid (based on the SPArse Matrix library), the libocn C library directly represents the network as a DAG. Each grid cell in the network has an associated outflow direction (given as an integer, 0-7, representing the 8 possible directions to neighboring cells) and a list of the directions of its neighbors (given as an 8-bit integer, where each bit indicates whether there is an edge connecting the cell to one of its 8 neighbors). libocn also implements functions to traverse and manipulate the network structure according to the simulated annealing algorithm described in the orginal paper.

# PyOCN
The PyOCN frontend is a Python package that provides a high-level interface to the libocn C library. PyOCN uses the NetworkX library to expose the network graph and provides additional functions to manipulate and analyze the graph, as well as export the graph to various formats, including GeoTIFF. The main class in PyOCN is the `OCN` class, which includes methods for constructing an initial configuration and running the optimization algorithm.

# Citing PyOCN
If you use PyOCN in your research, please cite this package and the original paper by Carraro et al.