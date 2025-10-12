
import matplotlib.pyplot as plt
import PyOCN as po
import numpy as np

ocn = po.OCN.from_net_type(
    net_type="E",
    dims=(16, 16),
    wrap=True,
    random_state=84712,
)
ocn.fit_custom_cooling(lambda t: np.ones_like(t)*1e-1, pbar=True, n_iterations=16*16*100, max_iterations_per_loop=1)
ocn.fit(max_iterations_per_loop=1)



# print(ocn.history.shape)
# print(np.max(np.diff(ocn.history[:, 1])))
# print(np.quantile(np.diff(ocn.history[:, 1]), 0.999))

plt.plot(ocn.history[:, 0], ocn.history[:, 1])
plt.plot(ocn.history[:, 0], ocn.history[:, 2])
# plt.xscale("log")
# plt.yscale("log")
plt.show()