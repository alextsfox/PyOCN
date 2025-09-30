#ifndef PYINIT_C
#define PYINIT_C

/*
setuptools.Extension expects the compiled library to export an initialization function
named PyINIT_<modulename>() so that the interpreter knows how to import it. Without it,
the import machinery cannot load the shared object as a python module. My raw libocn.so
compiled from gcc (as I use it to develop) is just a shared library, not a python module.
Adding this shim turns it into a module that can be imported.

In my case, I don't want to expose any python-callable functions from the C core. 
I just want the .so to exist inside the package to ctypes can load it. This shim
is the minimal valid python extension module definition. It defines a module, _libocn
that registers no methods. The result:
    * the .so passes python's extension loader checks
    * it is installed in the correct place (PyOCN/_libocn.cpyton-...so)
    * I can still load it with ctypes.CDLL bucause it's just a shared lib with extra python metadata.

So the shim is basically a tiny wrapper whose sole job is to satisfy Python’s extension 
import protocol and make your compiled C code compatible with the Extension build system.

If I wanted to, I could place my to_digraph constructor function here and expose it to python
using PyObjects, PyCapsule, whatever. This can squeeze out more performance and tighter integration.
If I take that approach, then libocn becomes part of the python runtime itself, not just a foreign 
binary playing dress-up as a python module (like I have now).
*/
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

#endif // PYINIT_C
