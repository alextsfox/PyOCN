from itertools import product

from tqdm import trange, tqdm
import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

import multiprocessing

import PyOCN


m, n = 50, 50
n_iterations = m*n*100
constant_phase = 0.1
cooling_rate = 0.2

def create_and_run(gamma):
    ocn = PyOCN.OCN.from_net_type(
        "H", 
        dims=(m, n), 
        random_state=8470,
        gamma=gamma,
        verbosity=0
    )
    energy = ocn.fit(energy_reports=1000, constant_phase=constant_phase, cooling_rate=cooling_rate, n_iterations=n_iterations, pbar=False)
    return ocn, energy

if __name__ == '__main__':
    gamma_range = [0.1, 0.5, 0.9]
    with multiprocessing.Pool(len(gamma_range)) as pool:
        results = list(pool.map(create_and_run, gamma_range))

    fig, axs = plt.subplots(2, 3, figsize=(15/2, 10/2), height_ratios=[0.5, 1])

    # high gamma
    ocn, energy = results[0]
    PyOCN.plot_ocn_energy_raster(ocn=ocn, norm=mpl.colors.LogNorm(), ax=axs[1,0])
    axs[0, 0].plot(np.linspace(0, n_iterations, len(energy)), energy/energy.max())


    # middling gamma
    ocn, energy = results[1]
    PyOCN.plot_ocn_energy_raster(ocn=ocn, norm=mpl.colors.LogNorm(), ax=axs[1,1])
    axs[0, 1].plot(np.linspace(0, n_iterations, len(energy)), energy/energy.max())

    # low gamma
    ocn, energy = results[2]
    PyOCN.plot_ocn_energy_raster(ocn=ocn, norm=mpl.colors.LogNorm(), ax=axs[1,2])
    axs[0, 2].plot(np.linspace(0, n_iterations, len(energy)), energy/energy.max())

    fig.tight_layout(h_pad=0.1)

    for ax in axs[1]: ax.set_axis_off()


    for ax, gamma in zip(axs[0], gamma_range):
        ax.set_ylim(0, 1)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_xlabel("Iteration")
        ax.set_ylabel("Normalized Energy")

        xticks = ax.get_xticks()
        ax.set_xticks(xticks, labels=[f"{int(x/1000)}k" for x in xticks])
        ax.set_xlim(0, n_iterations)
        ax.set_title(f"Gamma = {gamma:.2f}")

    plt.show()

