from itertools import product

from tqdm import trange, tqdm
import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

import PyOCN

ocn = PyOCN.OCN.from_net_type(
    "H", 
    dims=(90, 90), 
    random_state=8473
)

ocn.fit(max_iterations_per_loop=2_000, pbar=False)

PyOCN.plot_ocn_energy_raster(ocn=ocn, norm=mpl.colors.LogNorm(vmin=1, vmax=ocn.energy), cmap="terrain_r")
# # PyOCN.plot_ocn_as_dag(ocn, attribute='energy', with_labels=False, node_size=10)
plt.show()

