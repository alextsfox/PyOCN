from timeit import timeit

import PyOCN as po

ocn = po.OCN.from_net_type(
    net_type="V",
    dims=(74, 74),
    random_state=8472,
)
print(timeit(ocn.fit, number=10)/10)