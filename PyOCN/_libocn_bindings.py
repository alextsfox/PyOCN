"""
_StreamGraphC.py

Bindings for the StreamGraph C library, providing ctypes-based access to C data structures and functions
for stream graph manipulation and analysis.

Classes:
    _CartPairC    -- ctypes Structure mapping the C CartPairC struct
    _VertexC      -- ctypes Structure mapping the C Vertex struct
    _StreamGraphC -- ctypes Structure mapping the C StreamGraph struct


Usage:
    Loads the shared library 'libocn' and sets up function signatures for use in Python.

Author: Alexander S Fox
Copyright: (c) 2025 Alexander S Fox. All rights reserved.

This file is part of the OCN project.
"""

import numpy as np

from ctypes import (
    CDLL,
    Structure,
    c_uint16,
    c_uint32,
    c_double,
    c_uint8,
    c_bool,
    POINTER,
)
from pathlib import Path

CTYPE_TO_NP_DTYPE = {
    c_uint8: np.uint8,
    c_uint16: np.uint16,
    c_uint32: np.uint32,
    c_double: np.float64,
    c_bool: bool,
}

# TODO This will need to change probably depending on OS/architecture
_streamgraph_so_file = Path(__file__).parent.parent / "src" / "libocn.so"
libocn = CDLL(_streamgraph_so_file)


#############################
#   STATUS.H EQUIVALENTS    #
#############################
Status = c_uint8

SUCCESS = int(Status.in_dll(libocn, "SUCCESS").value)
OOB_ERROR = int(Status.in_dll(libocn, "OOB_ERROR").value)
NULL_POINTER_ERROR = int(Status.in_dll(libocn, "NULL_POINTER_ERROR").value)
SWAP_WARNING = int(Status.in_dll(libocn, "SWAP_WARNING").value)
MALFORMED_GRAPH_WARNING = int(Status.in_dll(libocn, "MALFORMED_GRAPH_WARNING").value)

STATUS_CODES = {
    SUCCESS: "SUCCESS",
    OOB_ERROR: "OOB_ERROR",
    NULL_POINTER_ERROR: "NULL_POINTER_ERROR",
    SWAP_WARNING: "SWAP_WARNING",
    MALFORMED_GRAPH_WARNING: "MALFORMED_GRAPH_WARNING",
}

#############################
# STREAMGRAPH.H EQUIVALENTS #
#############################
drainedarea_t = c_uint32
cartidx_t = c_uint16

class CartPair_C(Structure):
    _fields_ = [
        ("row", cartidx_t),
        ("col", cartidx_t),
    ]
    
linidx_t = c_uint32
localedges_t = c_uint8
clockhand_t = c_uint8
IS_ROOT = int(clockhand_t.in_dll(libocn, "IS_ROOT").value)


class Vertex_C(Structure):
    _fields_ = [
        ("drained_area", drainedarea_t),
        ("adown", linidx_t),
        ("edges", localedges_t),
        ("downstream", clockhand_t),
        ("visited", c_uint8),
    ]
vert_dtype = np.dtype([
    ("drained_area", CTYPE_TO_NP_DTYPE[drainedarea_t]),
    ("adown", CTYPE_TO_NP_DTYPE[linidx_t]),
    ("edges", CTYPE_TO_NP_DTYPE[localedges_t]),
    ("downstream", CTYPE_TO_NP_DTYPE[clockhand_t]),
    ("visited", CTYPE_TO_NP_DTYPE[c_uint8]),
])

class StreamGraph_C(Structure):
    _fields_ = [
        ("dims", CartPair_C),
        ("root", CartPair_C),
        ("energy", c_double),
        ("vertices", POINTER(Vertex_C)),  # Vertex*
    ]



# linidx_t sg_cart_to_lin(CartPairC coords, CartPairC dims);
libocn.sg_cart_to_lin.argtypes = [CartPair_C, CartPair_C]
libocn.sg_cart_to_lin.restype = linidx_t

# CartPairC sg_lin_to_cart(linidx_t a, CartPairC dims);
libocn.sg_lin_to_cart.argtypes = [linidx_t, CartPair_C]
libocn.sg_lin_to_cart.restype = CartPair_C

# Status sg_clockhand_to_lin_safe(linidx_t *a_down, linidx_t a, clockhand_t down, CartPairC dims);
libocn.sg_clockhand_to_lin_safe.argtypes = [POINTER(linidx_t), linidx_t, clockhand_t, CartPair_C]
libocn.sg_clockhand_to_lin_safe.restype = Status

# Status sg_get_cart_safe(Vertex *out, StreamGraph *G, CartPairC coords);
libocn.sg_get_cart_safe.argtypes = [POINTER(Vertex_C), POINTER(StreamGraph_C), CartPair_C]
libocn.sg_get_cart_safe.restype = Status

# Vertex sg_get_cart(StreamGraph *G, CartPairC coords);
libocn.sg_get_cart.argtypes = [POINTER(StreamGraph_C), CartPair_C]
libocn.sg_get_cart.restype = Vertex_C

# Status sg_set_cart_safe(StreamGraph *G, Vertex vert, CartPairC coords);
libocn.sg_set_cart_safe.argtypes = [POINTER(StreamGraph_C), Vertex_C, CartPair_C]
libocn.sg_set_cart_safe.restype = Status

# void sg_set_cart(StreamGraph *G, Vertex vert, CartPairC coords);
libocn.sg_set_cart.argtypes = [POINTER(StreamGraph_C), Vertex_C, CartPair_C]
libocn.sg_set_cart.restype = None

# Status sg_get_lin_safe(Vertex *out, StreamGraph *G, linidx_t a);
libocn.sg_get_lin_safe.argtypes = [POINTER(Vertex_C), POINTER(StreamGraph_C), linidx_t]
libocn.sg_get_lin_safe.restype = Status

# Vertex sg_get_lin(StreamGraph *G, linidx_t a);
libocn.sg_get_lin.argtypes = [POINTER(StreamGraph_C), linidx_t]
libocn.sg_get_lin.restype = Vertex_C

# Status sg_set_lin_safe(StreamGraph *G, Vertex vert, linidx_t a);
libocn.sg_set_lin_safe.argtypes = [POINTER(StreamGraph_C), Vertex_C, linidx_t]
libocn.sg_set_lin_safe.restype = Status

# void sg_set_lin(StreamGraph *G, Vertex vert, linidx_t a);
libocn.sg_set_lin.argtypes = [POINTER(StreamGraph_C), Vertex_C, linidx_t]
libocn.sg_set_lin.restype = None

# Status sg_create_empty_safe(StreamGraph *G, CartPairC root, CartPairC dims);
libocn.sg_create_empty_safe.argtypes = [POINTER(StreamGraph_C), CartPair_C, CartPair_C]
libocn.sg_create_empty_safe.restype = Status

# Status sg_destroy_safe(StreamGraph *G);
libocn.sg_destroy_safe.argtypes = [POINTER(StreamGraph_C)]
libocn.sg_destroy_safe.restype = Status

# Status sg_change_vertex_outflow(StreamGraph *G, linidx_t a, clockhand_t down_new);
libocn.sg_change_vertex_outflow.argtypes = [POINTER(StreamGraph_C), linidx_t, clockhand_t]
libocn.sg_change_vertex_outflow.restype = Status

# Status sg_flow_downstream_safe(StreamGraph *G, linidx_t a, uint8_t ncalls);
libocn.sg_flow_downstream_safe.argtypes = [POINTER(StreamGraph_C), linidx_t, c_uint8]
libocn.sg_flow_downstream_safe.restype = Status

# void sg_display(StreamGraph *G, bool use_utf8);
libocn.sg_display.argtypes = [POINTER(StreamGraph_C), c_bool]
libocn.sg_display.restype = None

# StreamGraph *sg_make_test_graph();
libocn.sg_make_test_graph.argtypes = []
libocn.sg_make_test_graph.restype = POINTER(StreamGraph_C)


##############################
#     OCN.H EQUIVALENTS      #
##############################

# Status ocn_update_energy(StreamGraph *G, drainedarea_t da_inc, linidx_t a, double gamma);
libocn.ocn_update_energy.argtypes = [POINTER(StreamGraph_C), drainedarea_t, linidx_t, c_double]
libocn.ocn_update_energy.restype = Status

# Status ocn_single_erosion_event(StreamGraph *G, double gamma, double temperature);
libocn.ocn_single_erosion_event.argtypes = [POINTER(StreamGraph_C), c_double, c_double]
libocn.ocn_single_erosion_event.restype = Status

# Status ocn_outer_ocn_loop(StreamGraph *G, uint32_t niterations, double gamma, double *annealing_schedule);
libocn.ocn_outer_ocn_loop.argtypes = [POINTER(StreamGraph_C), c_uint32, c_double, POINTER(c_double)]
libocn.ocn_outer_ocn_loop.restype = Status