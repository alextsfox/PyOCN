#include <Python.h>

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_libocn",
    "OCN C core (ctypes-loaded).",
    -1,
    NULL, NULL, NULL, NULL, NULL
};

PyMODINIT_FUNC PyInit__libocn(void) {
    return PyModule_Create(&moduledef);
}