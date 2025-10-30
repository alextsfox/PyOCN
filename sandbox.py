
import PyOCN as po
import ctypes
import matplotlib.pyplot as plt
from time import perf_counter as timer
from tqdm import trange
import numpy as np

ocn = po.OCN.from_net_type("H", dims=(10, 10))

pos = (7, 7)
print(list(ocn.predecessors(pos)))
print(list(ocn.successors(pos)))
print(po.utils.ancestors(ocn, pos))
print(po.utils.descendants(ocn, pos))

po.plotting.plot_ocn_as_dag(ocn)
plt.show()

