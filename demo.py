"""
Demonstration of generating Optimal Channel Networks (OCNs) using PyOCN with different gamma values.
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from time import perf_counter as timer

import multiprocessing

import PyOCN


m, n = 50, 50
n_iterations = m*n*40

def create_and_run(gamma):
    ocn = PyOCN.OCN.from_net_type(
        "H", 
        dims=(m, n), 
        random_state=8473,
        gamma=gamma,
        verbosity=0
    )
    energy = ocn.fit(
        constant_phase=0.01,
        cooling_rate=0.5, 
        energy_reports=1000, 
        n_iterations=n_iterations, 
        pbar=False
    )
    dag = ocn.to_digraph()
    # Multiprocessing can't pickle ocn objects, which contain pointers to C structs.
    # Instead, return the DAG and reconstruct the OCN later.
    # This would not be necessary if we were running single-threaded.
    return dag, energy  

if __name__ == '__main__':
    gamma_range = [0.1, 0.5, 0.9]
    with multiprocessing.Pool(len(gamma_range)) as pool:
        results = list(pool.map(create_and_run, gamma_range))
    print("Done generating OCNs, plotting...")

    fig, axs = plt.subplots(2, 3, figsize=(15/2, 8/2), height_ratios=[1, 0.25])

    # high gamma
    dag, energy = results[0]
    ocn = PyOCN.OCN.from_digraph(dag)  # Reconstruct OCN from DAG.
    PyOCN.plot_ocn_energy_raster(ocn=ocn, norm=mpl.colors.LogNorm(), ax=axs[0,0], cmap="cubehelix_r")
    axs[1, 0].plot(np.linspace(0, n_iterations, len(energy)), energy/energy.max())


    # middling gamma
    dag, energy = results[1]
    ocn = PyOCN.OCN.from_digraph(dag)
    PyOCN.plot_ocn_energy_raster(ocn=ocn, norm=mpl.colors.LogNorm(), ax=axs[0,1], cmap="cubehelix_r")
    axs[1, 1].plot(np.linspace(0, n_iterations, len(energy)), energy/energy.max())

    # low gamma
    dag, energy = results[2]
    ocn = PyOCN.OCN.from_digraph(dag)
    PyOCN.plot_ocn_energy_raster(ocn=ocn, norm=mpl.colors.LogNorm(), ax=axs[0,2], cmap="cubehelix_r")
    axs[1, 2].plot(np.linspace(0, n_iterations, len(energy)), energy/energy.max())

    fig.tight_layout(h_pad=0.5)

    for ax, gamma in zip(axs[0], gamma_range):
        ax.set_axis_off()
        ax.set_title(f"Gamma = {gamma:.2f}")


    for ax, gamma in zip(axs[1], gamma_range):
        ax.set_ylim(0, 1)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_xlabel("Iteration")
        ax.set_ylabel("Normalized Energy")

        xticks = ax.get_xticks()
        ax.set_xticks(xticks, labels=[f"{int(x/1000)}k" for x in xticks])
        ax.set_xlim(0, n_iterations)

    plt.show()

