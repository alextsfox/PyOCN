"""
Expose public exceptions & warnings
"""

import warnings

# Import status codes from _libocn_bindings
from ._libocn_bindings import (
    SUCCESS,
    OOB_ERROR,
    NULL_POINTER_ERROR,
    SWAP_WARNING,
    MALFORMED_GRAPH_WARNING,
)

STATUS_CODES = {
    SUCCESS: "SUCCESS",
    OOB_ERROR: "OOB_ERROR",
    NULL_POINTER_ERROR: "NULL_POINTER_ERROR",
    SWAP_WARNING: "SWAP_WARNING",
    MALFORMED_GRAPH_WARNING: "MALFORMED_GRAPH_WARNING",
}

class LibOCNError(Exception):
    """Base exception for libocn errors."""
    pass

class NullPointerError(LibOCNError):
    pass

class SwapWarning(RuntimeWarning):
    pass

class MalformedGraphWarning(RuntimeWarning):
    pass

STATUS_EXCEPTION_MAP = {
    OOB_ERROR: IndexError,
    NULL_POINTER_ERROR: NullPointerError,
}

STATUS_WARNING_MAP = {
    SWAP_WARNING: SwapWarning,
    MALFORMED_GRAPH_WARNING: MalformedGraphWarning,
}

def check_status(status):
    """
    Checks the status code returned by libocn functions.
    Raises an exception or issues a warning as appropriate.
    """
    if status == SUCCESS:
        return
    if status in STATUS_EXCEPTION_MAP:
        raise STATUS_EXCEPTION_MAP[status](f"libocn error: {STATUS_CODES.get(status, status)} ({status})")
    elif status in STATUS_WARNING_MAP:
        warnings.warn(f"libocn warning: {STATUS_CODES.get(status, status)} ({status})", STATUS_WARNING_MAP[status])
    else:
        raise LibOCNError(f"Unknown libocn status code: {status}")

__all__ = [
    "LibOCNError",
    "NullPointerError",
    "SwapWarning",
    "MalformedGraphWarning",
    "check_status",
]
