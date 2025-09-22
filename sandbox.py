import matplotlib.pyplot as plt
import networkx as nx

import PyOCN as ocn

sg = ocn.StreamGraph((3, 3), init_structure="test")

n_events = 4
fig, axs = plt.subplots(1, n_events, figsize=(12, 6))
sg.plot_streamgraph(ax=axs[0])
axs[0].set_title("Initial State, Energy {:.2f}".format(sg.energy))
for ax in axs[1:]:
    try:
        sg.single_erosion_event()
        label = f"SUCCESS {sg.energy:.2f}"
    except RuntimeError as e:
        label = f"FAILURE {sg.energy:.2f}"
        print(f"StreamGraph error during erosion event: {e}")
    sg.plot_streamgraph(ax=ax)
    ax.set_title(label)
plt.show()