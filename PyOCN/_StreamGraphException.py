import warnings
from ._libocn_bindings import libocn

class StreamGraphError(Exception): pass
class OutOfBoundsError(StreamGraphError): pass
class NullPointerError(StreamGraphError): pass

class StreamGraphWarning(RuntimeWarning): pass
class NoEdgeWarning(StreamGraphWarning): pass
class SwapWarning(StreamGraphWarning): pass
class MalformedGraphWarning(StreamGraphWarning): pass
class CycleWarning(StreamGraphWarning): pass

STATUS_EXCEPTION_MAP = {
    libocn.OOB_ERROR: OutOfBoundsError,
    libocn.NULL_POINTER_ERROR: NullPointerError,
}
STATUS_WARNING_MAP = {
    libocn.NO_EDGE_WARNING: NoEdgeWarning,
    libocn.SWAP_WARNING: SwapWarning,
    libocn.MALFORMED_GRAPH_WARNING: MalformedGraphWarning,
    libocn.CYCLE_WARNING: CycleWarning,
}

def check_status(status):
    if status == libocn.SUCCESS:
        return 
    if status in STATUS_EXCEPTION_MAP:
        raise STATUS_EXCEPTION_MAP[status](f"StreamGraph error code {status}")
    elif status in STATUS_WARNING_MAP:
        warnings.warn(f"StreamGraph warning code {status}", STATUS_WARNING_MAP[status])
    else:
        raise StreamGraphError(f"Unknown StreamGraph status code {status}")