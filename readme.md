This is a package to generate optimal channel networks, based on the algorithm described in Generation and application of river network analogues for use in ecology and evolution. Ecology and Evolution. doi:10.1002/ece3.6479 and on the OCNet R package (https://lucarraro.github.io/OCNet/).

# libocn
libocn is a C library that implements the core algorithms for generating and manipulating optimal channel networks (OCNs). It provides efficient data structures and functions to create, modify, and analyze OCNs.

Unlike the OCNet R package which uses an adjacency matrix implementation in R for manipulating the network (based on the SPARse Matrix library), the backend of PyOCN is the libocn C library, which uses
a custom data structure (FlowGrid) to represent the network. This structure represents the network graph directly, as a grid of cells with associated flow directions. Each cell has an associated outflow direction (given as an integer, 0-7, representing the 8 possible directions to neighboring cells) and a list of connected neighbors (given as an 8-bit integer, where each bit indicates whether there is a connection to the corresponding neighbor). This representation allows for efficient traversal and manipulation of the network, without the overhead of storing and manipulating large adjacency matrices. libocn also implements functions to traverse and manipulate the network structure according to the simulated annealing algorithm described in the orginal paper.

Optimizing a 512 x 512 grid to generate an OCN takes about 30 minutes on a laptop (MacBook Pro M1 Pro, 16GB RAM). Processing a 256 x 256 grid takes about 2 minutes.

# PyOCN
PyOCN is a Python package that provides a high-level interface to the libocn C library. It allows users to easily generate, manipulate, and analyze OCNs within Python. PyOCN uses the NetworkX library to expose the network graph and provides additional functions to export the graph in various formats, including as raster images.

<!-- # Citing PyOCN
If you use PyOCN in your research, please cite this software package as:

```
@software{fox_pyocn_2025,
  author = {Fox, Alexander S.},
  title = {{PyOCN}: A Python package for generating and manipulating optimal channel networks},
  url = {https://github.com/alexfox/pyocn},
  version = {0.1.0},
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

-->
