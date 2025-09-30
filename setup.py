from setuptools import setup, Extension, find_packages
import sys
from pathlib import Path


build_dir = Path("PyOCN")

c_src = build_dir / "c_src"
include_dirs = [str(c_src)]

sources = [
    c_src / "ocn.c",
    c_src / "flowgrid.c",
    c_src / "status.c",
    c_src / "rng.c",
    c_src / "pyinit.c",
]
sources = [str(s) for s in sources]

if sys.platform.startswith(("linux", "darwin")):
    extra_compile_args = ["-O3", "-flto", "-fPIC", "-std=c99", "-Wall", "-pedantic", "-march=native"]
    extra_link_args = ["-O3", "-flto"]
    # Link libm for pow() on Unix
    libraries = ["m"]
elif sys.platform.startswith("win"):
    extra_compile_args = ["/O2"]
    extra_link_args = []
    libraries = []
else:
    extra_compile_args = []
    extra_link_args = []
    libraries = []

# Define libocn as an extension module
# makes it possible to do `from PyOCN import libocn` if
# we ever decide to push expose the C API directly to the
# python runtime, instead of treating it as a foreign library.
ext = Extension(
    name="PyOCN._libocn",
    sources=sources,
    include_dirs=include_dirs,
    extra_compile_args=extra_compile_args,
    extra_link_args=extra_link_args,
    libraries=libraries,
)

setup(
    packages=find_packages(include=["PyOCN", "PyOCN.*"]),
    ext_modules=[ext],
    package_data={"PyOCN": []},  # this is where non-code or data files might go.
)