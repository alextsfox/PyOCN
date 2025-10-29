
from time import perf_counter as timer
import matplotlib.pyplot as plt
import PyOCN as po

rng = 8472
times = []
for _ in range(5):
    ocn = po.OCN.from_net_type(
        net_type="I",
        dims=(100, 100),
        wrap=True,
        random_state=rng,
    )
    timer_start = timer()
    ocn.fit(max_iterations_per_loop=10_000, pbar=True)
    timer_end = timer()
    rng += 1
    times.append(timer_end - timer_start)
print(f"Average time over 5 runs: {sum(times)/len(times):.2f} seconds")
print(f"Min time over 5 runs: {min(times):.2f} seconds")
print(f"Max time over 5 runs: {max(times):.2f} seconds")
print(f"Mean time over 5 runs: {sum(times)/len(times):.2f} seconds")
print(f"Final energy: {ocn.energy:.6f}")

plt.plot(ocn.history[:, 0], ocn.history[:, 1])
plt.plot(ocn.history[:, 0], ocn.history[:, 2])
plt.xscale("log")
plt.yscale("log")
plt.show()