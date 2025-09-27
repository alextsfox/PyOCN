# PyOCN
(c) 2024 Alexander S. Fox

This is a package to generate optimal channel networks (OCNs), based on the algorithm described in Carraro et al. (2020). *Generation and application of river network analogues for use in ecology and evolution. Ecology and Evolution.* doi:10.1002/ece3.6479 and mirrors much of the functionality of the OCNet R package (https://lucarraro.github.io/OCNet/).

This is a work in progress. I plan to release this as a package on PyPI in the near future. For now, you can clone the repository and build the package from source.

To compile libocn with `gcc` on MacOS or Linux, run the following commands from the root directory:

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

An initial stream network is generated as a directed acyclic graph (DAG) that is a spanning tree over a 2d grid of cells. Each cell in the grid (except the root) has a single outgoing edge that connects to one of its 8 neighboring cells. The OCN algorithm then iteratively modifies the DAG by randomly selecting a cell and changing its outflow to point to a different neighbor, ensuring that the spanning tree structure is maintained (*i.e. no cycles are introduced and each cell has a single outflow, except for the root node*). The total energy before and after the change is calculated based on the cumulative drained area:

$$
E[n] = \sum_{i} A[n]_i^\gamma
$$

Where $A[n]_i$ is the cumulative drained area at node $i$ at iteration $n$, and $\gamma$ is typically 0.5. This change is accepted or rejected based on simulated annealing, which is related to the Metropolis-Hastings algorithm. The probability of accepting a proposed change to the network is given by:

$$
P(accept) = \min(1, e^{-(E[n] - E[n-1]) / T[n]})
$$

Where $T[n]$ is the "temperature" of the network at iteration $n$. Initially, the temperature is set to a high value ($\sim E[0]$) to encourage exploration of the solution space. The temperature then exponentially decays over time, "annealing" the network as it settles into a low-energy configuration. Note that a proposed change is always accepted if it results in a lower energy state ($E[n] < E[n-1]$).

# libocn
The backend of PyOCN is the libocn C library. libocn implements the core algorithms for generating and manipulating OCNs. Unlike the OCNet R package which uses an adjacency matrix implementation in R for manipulating the network (based on the SPArse Matrix library), the libocn C library directly uses a directed acyclic graph (DAG) to represent the network. This structure represents the network as a grid of cells with associated flow directions. Each cell has an associated outflow direction (given as an integer, 0-7, representing the 8 possible directions to neighboring cells) and a list of the directions of its neighbors (given as an 8-bit integer, where each bit indicates whether there is a connection to the corresponding neighbor). This representation allows for efficient traversal and manipulation of the network, without the overhead of storing and manipulating large adjacency matrices. libocn also implements functions to traverse and manipulate the network structure according to the simulated annealing algorithm described in the orginal paper.

Optimizing a 512 x 512 grid to generate an OCN takes about 30 minutes on a laptop (MacBook Pro M1 Pro, 16GB RAM). Processing a 256 x 256 grid takes about 2 minutes.

# PyOCN
The PyOCN frontend is a Python package that provides a high-level interface to the libocn C library. It allows users to easily generate, manipulate, and analyze OCNs. PyOCN uses the NetworkX library to expose the network graph and provides additional functions to export the graph to various formats, including as raster images. (TODO: implement export to gtif, rasterio)

# Citing PyOCN
If you use PyOCN in your research, please cite this software package as:

```
@software{fox_pyocn_2025,
  author = {Fox, Alexander S.},
  title = {{PyOCN and libocn}},
  url = {https://github.com/alexfox/pyocn},
  date = {2025-10-01},
} 
```

in addition to the original paper:

```
@article{carraro_generation_2020,
	title = {Generation and application of river network analogues for use in ecology and evolution},
	volume = {10},
	doi = {10.1002/ece3.6479},
	number = {14},
	journal = {Ecology and Evolution},
	author = {Carraro, Luca and Bertuzzo, Enrico and Fronhofer, Emanuel A. and Furrer, Reinhard and Gounand, Isabelle and Rinaldo, Andrea and Altermatt, Florian},
	year = {2020},
	pages = {7537--7550},
}
```
