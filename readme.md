# OCN
Program to generate optimal channel networks.

https://onlinelibrary-wiley-com.uwyo.idm.oclc.org/doi/full/10.1002/ece3.6479

Let us consider a regular lattice made up of $N$ cells, where each cell represents the generic node $i$ of the network. Each node $i$ is connected via a link to one of its nearest neighbors. The energy dissipation across the $i$th network link is proportional to $Q_i \Delta h_i$, where $Q_i$ is the landscape-forming discharge in the link, and $\Delta h_i = s_i L_i$ the corresponding elevation drop, with $s_i$ identifying slope and $L_i$ link length. By assuming $Q_i \sim A_i$, where $A_i$ is the area drained by link $i$, and the slope–area relationship $s_i \sim A^{\gamma - 1}$, the functional representing total energy expenditure across a landscape formed by $N$ cells reads.

$$
H = \sum_{i}^N A_i^\gamma \tag{1}
$$

Note that link lengths do not appear in the above formula, as they can be considered constant with no loss of generality. The OCN configuration is defined by an adjacency matrix $W$ whose corresponding set of drainage areas $A=[A_1,\dots,A_n]$ yields a local, dynamically accessible minimum of Equation (1). Note that the correspondence between $A$ and the adjacency matrix $W$ of a tree is subsumed by the relationship $I_N-W^TA=1$, where $I_N$ is the identity matrix of order $N$, and 1 a $N$-dimensional vector of ones.

Minimization of Equation (1) is operated by means of a simulated annealing technique: starting from a feasible initial flow configuration (i.e., a spanning tree), a link at a time is rewired to one of its nearest neighbors; if the obtained configuration is a spanning tree, $H$ is computed; the new configuration is accepted if it lowers total energy expenditure; if this is not the case, the new configuration can still be accepted with a probability controlled by the cooling schedule of the simulated annealing algorithm. Such myopic search, which only explores close configurations, actually mimics the type of optimization that nature performs, at least in fluvial landscapes. Notably, restricting the search of a network yielding a minimum of Equation (1) to spanning, loopless configurations entails no approximation, because every spanning tree is a local minimum of total energy dissipation. The shape of the so-obtained OCN retains the heritage of the initial flow configuration, although the extent to which this occurs is partly controlled by the cooling schedule adopted. This underpins the concept of feasible optimality, that is, the search for dynamically accessible configurations.

Coding challanges:
* Implement a 2D Array structre in C with basic operations (addition, subtraction, multiplication, division, negation, element access, setting, slicing, copying) (DONE)
* Implement a spanning tree structure in C with basic operations (adding/removing edges, checking for cycles, finding paths, swapping/permuting edges)
* Implement the OCN algorithm using simulated annealing to minimize the energy functional H
* Implement input/output functions to read/write arrays and trees from/to files