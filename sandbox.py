import warnings
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

import PyOCN


def to_c_code(ocn:PyOCN.OCN, label):
    # StreamGraph sg;
    print(f"// {label}")
    print("StreamGraph sg;")
    print("CartPair dims = {", ocn.sg.dims[0], ", ", ocn.sg.dims[1], "};", sep="")
    print("CartPair root = {", ocn.sg.root[0], ", ", ocn.sg.root[1], "};", sep="")
    print("sg_create_empty_safe(&sg, root, dims);")
    print("Vertex vertices[", ocn.sg.size, "] = {", sep="")
    for i, v in enumerate(ocn.sg):
        print("    {", end="")
        print(".drained_area =", v.drained_area, end=", ")
        print(".adown =", v.adown, end=", ")
        print(".edges =", v.edges, end=", ")
        print(".downstream =", v.downstream, end=", ")
        print(".visited =", v.visited, end="")
        if i < ocn.sg.size - 1: print("},")
        else: print("}")
    print("};")

    print("sg.vertices = vertices;")
    print("sg.energy =", ocn.energy, ";")


################################################
# basic 4x4 graph, no issues
label = "basic 4x4"
G = nx.DiGraph()

G.add_node(0, pos=(0, 0), drained_area=1)
G.add_node(1, pos=(0, 1), drained_area=1)
G.add_node(2, pos=(0, 2), drained_area=1)
G.add_node(3, pos=(0, 3), drained_area=1)
G.add_node(4, pos=(1, 0), drained_area=2)
G.add_node(5, pos=(1, 1), drained_area=2)
G.add_node(6, pos=(1, 2), drained_area=2)
G.add_node(7, pos=(1, 3), drained_area=2)
G.add_node(8, pos=(2, 0), drained_area=3)
G.add_node(9, pos=(2, 1), drained_area=3)
G.add_node(10, pos=(2, 2), drained_area=3)
G.add_node(11, pos=(2, 3), drained_area=3)
G.add_node(12, pos=(3, 0), drained_area=4)
G.add_node(13, pos=(3, 1), drained_area=8)
G.add_node(14, pos=(3, 2), drained_area=12)
G.add_node(15, pos=(3, 3), drained_area=16)

G.add_edge(0, 4)
G.add_edge(1, 5)
G.add_edge(2, 6)
G.add_edge(3, 7)
G.add_edge(4, 8)
G.add_edge(5, 9)
G.add_edge(6, 10)
G.add_edge(7, 11)
G.add_edge(8, 12)
G.add_edge(9, 13)
G.add_edge(10, 14)
G.add_edge(11, 15)
G.add_edge(12, 13)
G.add_edge(13, 14)
G.add_edge(14, 15)

ocn = PyOCN.OCN(init_structure=G, gamma=1, random_state=8472)

to_c_code(ocn, label)

# ############################################
# # contains an "X"
# label = "X shape, 5->10 and 6->9"
# G = nx.DiGraph()

# G.add_node(0, pos=(0, 0), drained_area=1)
# G.add_node(1, pos=(0, 1), drained_area=1)
# G.add_node(2, pos=(0, 2), drained_area=1)
# G.add_node(3, pos=(0, 3), drained_area=1)
# G.add_node(4, pos=(1, 0), drained_area=2)
# G.add_node(5, pos=(1, 1), drained_area=2)
# G.add_node(6, pos=(1, 2), drained_area=2)
# G.add_node(7, pos=(1, 3), drained_area=2)
# G.add_node(8, pos=(2, 0), drained_area=3)
# G.add_node(9, pos=(2, 1), drained_area=3)
# G.add_node(10, pos=(2, 2), drained_area=3)
# G.add_node(11, pos=(2, 3), drained_area=3)
# G.add_node(12, pos=(3, 0), drained_area=4)
# G.add_node(13, pos=(3, 1), drained_area=8)
# G.add_node(14, pos=(3, 2), drained_area=12)
# G.add_node(15, pos=(3, 3), drained_area=16)

# G.add_edge(0, 4)
# G.add_edge(1, 5)
# G.add_edge(2, 6)
# G.add_edge(3, 7)
# G.add_edge(4, 8)
# G.add_edge(5, 10)
# G.add_edge(6, 9)
# G.add_edge(7, 11)
# G.add_edge(8, 12)
# G.add_edge(9, 13)
# G.add_edge(10, 14)
# G.add_edge(11, 15)
# G.add_edge(12, 13)
# G.add_edge(13, 14)
# G.add_edge(14, 15)

# ocn = PyOCN.OCN(init_structure=G, gamma=1, random_state=8472)
# print(ocn.sg._vertices.flatten()[9])

# to_c_code(ocn, label)


# ############################################
# # contains a cycle
# G = nx.DiGraph()

# G.add_node(0, pos=(0, 0), drained_area=1)
# G.add_node(1, pos=(0, 1), drained_area=1)
# G.add_node(2, pos=(0, 2), drained_area=1)
# G.add_node(3, pos=(0, 3), drained_area=1)
# G.add_node(4, pos=(1, 0), drained_area=2)
# G.add_node(5, pos=(1, 1), drained_area=2)
# G.add_node(6, pos=(1, 2), drained_area=2)
# G.add_node(7, pos=(1, 3), drained_area=2)
# G.add_node(8, pos=(2, 0), drained_area=3)
# G.add_node(9, pos=(2, 1), drained_area=3)
# G.add_node(10, pos=(2, 2), drained_area=3)
# G.add_node(11, pos=(2, 3), drained_area=3)
# G.add_node(12, pos=(3, 0), drained_area=4)
# G.add_node(13, pos=(3, 1), drained_area=8)
# G.add_node(14, pos=(3, 2), drained_area=12)
# G.add_node(15, pos=(3, 3), drained_area=16)

# G.add_edge(0, 4)
# G.add_edge(1, 5)
# G.add_edge(2, 6)
# G.add_edge(3, 7)
# G.add_edge(4, 8)
# G.add_edge(5, 6)
# G.add_edge(6, 10)
# G.add_edge(7, 11)
# G.add_edge(8, 12)
# G.add_edge(9, 5)
# G.add_edge(10, 9)
# G.add_edge(11, 15)
# G.add_edge(12, 13)
# G.add_edge(13, 14)
# G.add_edge(14, 15)

# PyOCN.plot_positional_digraph(G)

"""
// contains a cycle
StreamGraph sg;
CartPair dims = {4, 4};
CartPair root = {3, 3};
sg_create_empty_safe(&sg, root, dims);
Vertex vertices[16] = {
    {.drained_area = 1, .adown = 4, .edges = 16, .downstream = 4, .visited = 0},
    {.drained_area = 1, .adown = 5, .edges = 16, .downstream = 4, .visited = 0},
    {.drained_area = 1, .adown = 6, .edges = 16, .downstream = 4, .visited = 0},
    {.drained_area = 1, .adown = 7, .edges = 16, .downstream = 4, .visited = 0},
    {.drained_area = 2, .adown = 8, .edges = 17, .downstream = 4, .visited = 0},
    {.drained_area = 2, .adown = 6, .edges = 21, .downstream = 4, .visited = 0},  // 5 goes to 6
    {.drained_area = 2, .adown = 10, .edges = 81, .downstream = 4, .visited = 0},  // 6 goes to 10
    {.drained_area = 2, .adown = 11, .edges = 17, .downstream = 4, .visited = 0},
    {.drained_area = 3, .adown = 12, .edges = 17, .downstream = 4, .visited = 0},
    {.drained_area = 3, .adown = 5, .edges = 5, .downstream = 4, .visited = 0},  // 9 goes to 5
    {.drained_area = 3, .adown = 9, .edges = 65, .downstream = 4, .visited = 0},  // 10 goes to 9
    {.drained_area = 3, .adown = 15, .edges = 17, .downstream = 4, .visited = 0},
    {.drained_area = 4, .adown = 13, .edges = 5, .downstream = 2, .visited = 0},
    {.drained_area = 8, .adown = 14, .edges = 69, .downstream = 2, .visited = 0},
    {.drained_area = 12, .adown = 15, .edges = 69, .downstream = 2, .visited = 0},
    {.drained_area = 16, .adown = 0, .edges = 65, .downstream = 255, .visited = 0}
};
sg.vertices = vertices;
sg.energy = 64.0 ;
"""

# ################################
# # contains staring contest
# """
# // staring contest
# StreamGraph sg;
# CartPair dims = {4, 4};
# CartPair root = {3, 3};
# sg_create_empty_safe(&sg, root, dims);
# Vertex vertices[16] = {
#     {.drained_area = 1, .adown = 4, .edges = 16, .downstream = 4, .visited = 0},
#     {.drained_area = 1, .adown = 5, .edges = 16, .downstream = 4, .visited = 0},
#     {.drained_area = 1, .adown = 6, .edges = 16, .downstream = 4, .visited = 0},
#     {.drained_area = 1, .adown = 7, .edges = 16, .downstream = 4, .visited = 0},
#     {.drained_area = 2, .adown = 8, .edges = 16, .downstream = 4, .visited = 0},
#     {.drained_area = 2, .adown = 9, .edges = 16, .downstream = 4, .visited = 0},  // 5 looks at 9
#     {.drained_area = 2, .adown = 10, .edges = 16, .downstream = 4, .visited = 0},
#     {.drained_area = 2, .adown = 11, .edges = 16, .downstream = 4, .visited = 0},
#     {.drained_area = 3, .adown = 12, .edges = 16, .downstream = 4, .visited = 0},
#     {.drained_area = 3, .adown = 5, .edges = 16, .downstream = 4, .visited = 0},  // 9 looks at 5
#     {.drained_area = 3, .adown = 14, .edges = 16, .downstream = 4, .visited = 0},
#     {.drained_area = 3, .adown = 15, .edges = 16, .downstream = 4, .visited = 0},
#     {.drained_area = 4, .adown = 13, .edges = 4, .downstream = 2, .visited = 0},
#     {.drained_area = 8, .adown = 14, .edges = 4, .downstream = 2, .visited = 0},
#     {.drained_area = 12, .adown = 15, .edges = 4, .downstream = 2, .visited = 0},
#     {.drained_area = 16, .adown = 0, .edges = 0, .downstream = 255, .visited = 0}
# };
# sg.vertices = vertices;
# sg.energy = 64.0 ;
# """

# plt.show()
