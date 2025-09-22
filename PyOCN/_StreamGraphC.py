"""
_StreamGraphC.py

Bindings for the StreamGraph C library, providing ctypes-based access to C data structures and functions
for stream graph manipulation and analysis.

Classes:
    _VertexC      -- ctypes Structure mapping the C Vertex struct
    _StreamGraphC -- ctypes Structure mapping the C StreamGraph struct

Usage:
    Loads the shared library 'streamgraph.so' and sets up function signatures for use in Python.

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
    POINTER,
)
from pathlib import Path

class VertexC(Structure):
    _fields_ = [
        ("drained_area", c_uint32),  # drainedarea_t
        ("adown", c_uint32),         # linidx_t
        ("edges", c_uint8),          # localedges_t
        ("downstream", c_uint8),     # clockhand_t
        ("visited", c_uint8),        # uint8_t
    ]
vert_dtype = np.dtype([
    ("drained_area", np.uint32),
    ("adown", np.uint32),
    ("edges", np.uint8),
    ("downstream", np.uint8),
    ("visited", np.uint8),
])

class StreamGraphC(Structure):
    _fields_ = [
        ("m", c_uint16),             # cartidx_t
        ("n", c_uint16),             # cartidx_t
        ("i_root", c_uint16),        # cartidx_t
        ("j_root", c_uint16),        # cartidx_t
        ("energy", c_double),
        ("vertices", POINTER(VertexC)),  # Vertex*
    ]

class cartesian_pair(Structure):
    _fields_ = [
        ("i", c_uint16),  # cartidx_t
        ("j", c_uint16),  # cartidx_t
    ]

_streamgraph_so_file = Path(__file__).parent.parent / "src" / "libocn.so"
libocn = CDLL(_streamgraph_so_file)

libocn.sg_make_test_graph.restype = StreamGraphC
libocn.sg_single_erosion_event.restype = c_uint8  # Status
libocn.sg_outer_ocn_loop.restype = c_uint8  # Status
libocn.sg_lin_to_cart.restype = cartesian_pair
libocn.sg_cart_to_lin.restype = c_uint32  # linidx_t

STATUS_CODES = {
    0: "SUCCESS",
    1: "OOB_ERROR",
    2: "NO_EDGE_WARNING",
    3: "NULL_POINTER_ERROR",
    4: "SWAP_WARNING",
    5: "MALFORMED_GRAPH_WARNING",
    6: "CYCLE_WARNING",
}

IS_ROOT = 255