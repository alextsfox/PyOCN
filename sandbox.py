
import PyOCN as po
import ctypes
import matplotlib.pyplot as plt
from time import perf_counter as timer
from tqdm import trange
import numpy as np

# def check_upstream_function_result(numpstream_iter, upstream_indices_iter, correct_answer, status):
#     # print("Iterative:")
#     # print("\tNumber of upstream vertices:", numpstream_iter)
#     # print("Correct answer:")
#     # print("\tNumber of upstream vertices:", len(correct_answer))
#     # print("Status code:", status)
#     success = (numpstream_iter == len(correct_answer)) and (upstream_indices_iter == correct_answer) and (status == 0)
#     # print(f"Success: {success}")
#     return success

# def run_upstream_functions(ocn, a):
#     m, n = ocn.dims
#     a = po._libocn_bindings.linidx_t(a)

#     c_graph = ocn._OCN__p_c_graph
#     upstream_indices = (po._libocn_bindings.linidx_t * (m * n))()
#     stack = (po._libocn_bindings.linidx_t * (m * n))()
#     nupstream_val = po._libocn_bindings.linidx_t(0)

#     t0 = timer()
#     status = po._libocn_bindings.libocn.fg_dfs_iterative(
#         upstream_indices,
#         ctypes.byref(nupstream_val),
#         stack,
#         c_graph,
#         a
#     )
#     t1 = timer()
#     po._libocn_bindings.libocn.ocn_compute_energy(c_graph, ocn.gamma)
#     t2 = timer()

#     numpstream_iter = nupstream_val.value
#     upstream_indices_iter = sorted(list(upstream_indices)[:nupstream_val.value])

#     return numpstream_iter, upstream_indices_iter, status, t1 - t0, t2 - t1

# ocn = po.OCN.from_net_type("I", dims=(10, 10))
# a = 85
# correct_answer = [80, 81, 82, 83, 84, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99]
# numpstream_iter, upstream_indices_iter, status, time_dfs, time_energy = run_upstream_functions(ocn, a)
# check_upstream_function_result(numpstream_iter, upstream_indices_iter, correct_answer, status)

# ocn = po.OCN.from_net_type("I", dims=(10, 10))
# a = 33
# correct_answer = [30, 31, 32]
# numpstream_iter, upstream_indices_iter, status, time_dfs, time_energy = run_upstream_functions(ocn, a)
# check_upstream_function_result(numpstream_iter, upstream_indices_iter, correct_answer, status)


ocn = po.OCN.from_net_type("E", dims=(64, 64), random_state=8471)
ocn.fit(pbar=True, max_iterations_per_loop=100, calculate_full_energy=True)
print("Final energy:", ocn.energy)

ocn = po.OCN.from_net_type("E", dims=(64, 64), random_state=8471)
ocn.fit(pbar=True, max_iterations_per_loop=100, calculate_full_energy=False)
print("Final energy:", ocn.energy)
# plt.hist(np.diff(ocn.history[:, 1]), bins=np.linspace(-1000, 0))

plt.show()