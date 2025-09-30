from setuptools import setup, Extension, find_packages
import sys
from pathlib import Path


sources = [
    str(Path("PyOCN")/"c_src"/"pymodule_shim.c"),
    str(Path("PyOCN")/"c_src"/"ocn.c"),
    str(Path("PyOCN")/"c_src"/"flowgrid.c"),
    str(Path("PyOCN")/"c_src"/"status.c"),
    str(Path("PyOCN")/"c_src"/"rng.c"),
]
include_dirs = [str(Path("PyOCN")/"c_src")]

extra_compile_args = []
extra_link_args = []

if sys.platform.startswith(("linux", "darwin")):
    extra_compile_args += ["-O3", "-flto", "-fPIC", "-std=c99", "-Wall", "-pedantic"]
    extra_link_args += ["-O3", "-flto"]
    # Link libm for pow() on Unix
    libraries = ["m"]
elif sys.platform.startswith("win"):
    # MSVC flags (optimize)
    extra_compile_args += ["/O2"]
    libraries = []  # msvcrt provides math
else:
    libraries = []

ext = Extension(
    name="PyOCN._libocn",                # Built extension inside the package
    sources=sources,
    include_dirs=include_dirs,
    extra_compile_args=extra_compile_args,
    extra_link_args=extra_link_args,
    libraries=libraries,
)

setup(
    packages=find_packages(include=["PyOCN", "PyOCN.*"]),
    ext_modules=[ext],
    package_data={"PyOCN": []},  # wheels will contain the built extension
)