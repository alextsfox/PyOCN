"""
libocn exception handling
"""
#TODO: I feel like this is extremely clunky and and doens't really give very helpful errors.

import warnings

# Import status codes from _libocn_bindings
from ._libocn_bindings import (
    SUCCESS,
    OOB_ERROR,
    EROSION_FAILURE,
    NULL_POINTER_ERROR,
    SWAP_WARNING,
    MALFORMED_GRAPH_WARNING,
)

STATUS_CODES = {
    SUCCESS: "SUCCESS",
    EROSION_FAILURE: "EROSION_FAILURE",
    OOB_ERROR: "OOB_ERROR",
    NULL_POINTER_ERROR: "NULL_POINTER_ERROR",
    SWAP_WARNING: "SWAP_WARNING",
    MALFORMED_GRAPH_WARNING: "MALFORMED_GRAPH_WARNING",
}

class LibOCNError(Exception):
    """Base exception for libocn errors."""
    pass

class NullPointerError(LibOCNError):
    default_message = (
        "\nTrying to access a null pointer. "
        "Something has gone terribly wrong inside PyOCN or libocn. "
        "Please open an issue at https://github.com/alextsfox/PyOCN/issues."
    )

    def __init__(self, custom_message=None):
        self.custom_message = custom_message

    def __str__(self):
        if self.custom_message:
            return f"{self.default_message}\nDetails: {self.custom_message}"
        return self.default_message

class SwapWarning(RuntimeWarning):
    pass

class MalformedGraphWarning(RuntimeWarning):
    default_message = (
        "\nYour flowgrid is invalid. "
        "This can be due to either cycles in the graph, "
        "the root node not having access to all nodes, "
        "incorrect dimensions (must have even rows and columns), "
        "or other structural issues."
    )

    def __init__(self, custom_message=None):
        self.custom_message = custom_message

    def __str__(self):
        if self.custom_message:
            return f"{self.default_message}\nDetails: {self.custom_message}"
        return self.default_message

STATUS_EXCEPTION_MAP = {
    OOB_ERROR: IndexError,
    NULL_POINTER_ERROR: NullPointerError,
}

STATUS_WARNING_MAP = {
    SWAP_WARNING: SwapWarning,
    MALFORMED_GRAPH_WARNING: MalformedGraphWarning,
}
STATUS_OKAY_MAP = {
    EROSION_FAILURE: None, 
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
    elif status in STATUS_OKAY_MAP:
        return
    else:
        raise LibOCNError(f"Unknown libocn status code: {status}")

__all__ = [
    "LibOCNError",
    "NullPointerError",
    "SwapWarning",
    "MalformedGraphWarning",
    "check_status",
]
