"""
Bindings for the libocn C library, 
providing ctypes-based access to C data structures and functions.
"""

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
import importlib
from pathlib import Path

def _load_libocn():
    # attempt the dev build first
    _dev_lib = Path(__file__).parent / "libocn_dev.so"
    if _dev_lib.exists():
        return CDLL(str(_dev_lib))
    
    # prefer the built extension inside the wheel: PyOCN/_libocn*.so|.pyd
    try:
        mod = importlib.import_module("._libocn", package=__package__)
        return CDLL(str(Path(mod.__file__)))
    except Exception:
        pass

    raise OSError(
        "Could not locate libocn: neither the built extension (PyOCN._libocn) " \
        "nor a local shared library was found. Reinstall the wheel or build " \
        "the C library manually with PyOCN/c_src/build.sh."
    )

libocn = _load_libocn()

#############################
#   STATUS.H EQUIVALENTS    #
#############################
Status = c_uint8

SUCCESS = int(Status.in_dll(libocn, "SUCCESS").value)
EROSION_FAILURE = int(Status.in_dll(libocn, "EROSION_FAILURE").value)
OOB_ERROR = int(Status.in_dll(libocn, "OOB_ERROR").value)
NULL_POINTER_ERROR = int(Status.in_dll(libocn, "NULL_POINTER_ERROR").value)
SWAP_WARNING = int(Status.in_dll(libocn, "SWAP_WARNING").value)
MALFORMED_GRAPH_WARNING = int(Status.in_dll(libocn, "MALFORMED_GRAPH_WARNING").value)

#############################
# STREAMGRAPH.H EQUIVALENTS #
#############################
drainedarea_t = c_double
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

class FlowGrid_C(Structure):
    _fields_ = [
        ("dims", CartPair_C),
        ("energy", c_double),
        ("resolution", c_double),
        ("nroots", c_uint16),
        ("vertices", POINTER(Vertex_C)),  # Vertex*
        ("wrap", c_bool),
    ]

# linidx_t fg_cart_to_lin_export(CartPairC coords, CartPairC dims);
libocn.fg_cart_to_lin_export.argtypes = [CartPair_C, CartPair_C]
libocn.fg_cart_to_lin_export.restype = linidx_t

# Status fg_get_lin_export(Vertex *out, FlowGrid *G, linidx_t a);
libocn.fg_get_lin_export.argtypes = [POINTER(Vertex_C), POINTER(FlowGrid_C), linidx_t]
libocn.fg_get_lin_export.restype = Status

# Status fg_set_lin_export(FlowGrid *G, Vertex vert, linidx_t a);
libocn.fg_set_lin_export.argtypes = [POINTER(FlowGrid_C), Vertex_C, linidx_t]
libocn.fg_set_lin_export.restype = Status

# Status fg_create_empty(CartPairC dims);
libocn.fg_create_empty.argtypes = [CartPair_C]
libocn.fg_create_empty.restype = POINTER(FlowGrid_C)

# FlowGrid *fg_copy(FlowGrid *G);
libocn.fg_copy.argtypes = [POINTER(FlowGrid_C)]
libocn.fg_copy.restype = POINTER(FlowGrid_C)

# Status fg_destroy(FlowGrid *G);
libocn.fg_destroy.argtypes = [POINTER(FlowGrid_C)]
libocn.fg_destroy.restype = Status

##############################
#     OCN.H EQUIVALENTS      #
##############################

# double ocn_compute_energy(FlowGrid *G, double gamma);
libocn.ocn_compute_energy.argtypes = [POINTER(FlowGrid_C), c_double]
libocn.ocn_compute_energy.restype = c_double

# Status ocn_single_erosion_event(FlowGrid *G, double gamma, double temperature);
libocn.ocn_single_erosion_event.argtypes = [POINTER(FlowGrid_C), c_double, c_double]
libocn.ocn_single_erosion_event.restype = Status

# Status ocn_outer_ocn_loop(FlowGrid *G, uint32_t niterations, double gamma, double *annealing_schedule);
libocn.ocn_outer_ocn_loop.argtypes = [POINTER(FlowGrid_C), c_uint32, c_double, POINTER(c_double)]
libocn.ocn_outer_ocn_loop.restype = Status

########################
#  RNG.H EQUIVALENTS   #
########################

# void rng_seed(unsigned int seed);
libocn.rng_seed.argtypes = [c_uint32]
libocn.rng_seed.restype = None

# void rng_seed_random();
libocn.rng_seed_random.argtypes = []
libocn.rng_seed_random.restype = None

__all__ = [
    "SUCCESS",
    "EROSION_FAILURE",
    "OOB_ERROR",
    "NULL_POINTER_ERROR",
    "SWAP_WARNING",
    "MALFORMED_GRAPH_WARNING",
    "libocn",
    "CartPair_C",
    "Vertex_C",
    "FlowGrid_C",
    "IS_ROOT",
]