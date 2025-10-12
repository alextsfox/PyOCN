
import matplotlib.pyplot as plt
import PyOCN as po

ocn = po.OCN.from_net_type(
    net_type="V",
    dims=(100, 100),
    wrap=True,
    random_state=8472,
)
ocn.fit(pbar=True, n_iterations=100**2*100, constant_phase=0.25, cooling_rate=1)

plt.plot(ocn.history[:, 0], ocn.history[:, 1])
plt.plot(ocn.history[:, 0], ocn.history[:, 2])
plt.show()