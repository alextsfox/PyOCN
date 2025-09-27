from itertools import product

from tqdm import trange, tqdm
import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

import PyOCN

ocn = PyOCN.OCN.from_net_type(
    "H", 
    dims=(256, 256), 
    random_state=8473,
    verbosity=2
)

energy = ocn.fit(report_energy_interval=3000)
plt.plot(energy)

# PyOCN.plot_ocn_energy_raster(ocn=ocn, norm=mpl.colors.LogNorm(vmin=1, vmax=ocn.energy), cmap="terrain_r")
# # PyOCN.plot_ocn_as_dag(ocn, attribute='energy', with_labels=False, node_size=10)
plt.show()

