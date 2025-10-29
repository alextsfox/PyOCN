
from time import perf_counter as timer
from concurrent.futures import ThreadPoolExecutor

import matplotlib.pyplot as plt
import numpy as np
import PyOCN as po


rng = 8472
times = []

gammas = np.linspace(0, 1, 51)
ocns = [po.OCN.from_net_type(
    net_type="I",
    gamma=gamma,
    dims=(100, 100),
    wrap=True,
    random_state=rng,
) for gamma in gammas]

fit_kwargs = {
    "n_iterations": 64*64*40,
    "pbar": False,
    "max_iterations_per_loop": 2_000,
}
results, fitted_ocns = po.utils.multi_fit(ocns, len(gammas), pbar=True, fit_kwargs=fit_kwargs)

for fitted_ocn in fitted_ocns:
    plt.plot(fitted_ocn.history[:, 0], (fitted_ocn.history[:, 1] - fitted_ocn.history[-1, 1]) / (fitted_ocn.history[0, 1] - fitted_ocn.history[-1, 1]), alpha=0.5, color="C0")
# print(f"Average time over 5 runs: {sum(times)/len(times):.2f} seconds")
# print(f"Min time over 5 runs: {min(times):.2f} seconds")
# print(f"Max time over 5 runs: {max(times):.2f} seconds")
# print(f"Mean time over 5 runs: {sum(times)/len(times):.2f} seconds")
# print(f"Final energy: {ocn.energy:.6f}")

plt.xscale("log")
plt.yscale("log")
plt.show()