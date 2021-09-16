from definitions import PY_SRC_DIR, C_SRC_DIR, C_INC_DIR
from distutils.core import setup, Extension
from sysconfig import get_paths
from pathlib import Path
import numpy as np
import glob
import os


def get_cpp_sources():
    srcs = [file for file in sorted(glob.glob(f"{C_SRC_DIR}/*.cpp"))]
    to_remove = ["write_ics.cpp", "main.cpp"]
    for src_file in srcs:
        for r in to_remove:
            if r in src_file:
                srcs.remove(src_file)

    return srcs


def get_py_sources():
    srcs = [Path(file).stem for file in sorted(glob.glob(f"{PY_SRC_DIR}/*.py"))]
    while "__init__" in srcs:
        srcs.remove("__init__")

    return srcs


# Getting Python and numpy src dirs
info = get_paths()
python_incl = info["include"]
python_lib = info["stdlib"]
numpy_incl = os.path.join(np.get_include(), "numpy")

# NCI
# install_py_modules = True
# armadillo_incl = "/home/149/ub0692/include/"
# armadillo_lib = "/home/149/ub0692/lib64"
# hdf5_incl = "/apps/hdf5/1.10.5/include"

# Home
install_py_modules = False
hdf5_incl = "/usr/include/hdf5/serial/"
hdf5_lib_dir = "/usr/lib/x86_64-linux-gnu/hdf5/serial/"
armadillo_incl = "/"
armadillo_lib = "/"


# Source files and definitions for this thing
include_dirs = [python_incl, numpy_incl, armadillo_incl, hdf5_incl, C_INC_DIR]
library_dirs = [python_lib, armadillo_lib, hdf5_lib_dir]
libraries = ["gsl", "gslcblas", "m", "armadillo", "hdf5"]
define_macros = [("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION"), ('H5_BUILT_AS_DYNAMIC_LIB', True)]
sources = get_cpp_sources()

if install_py_modules:
    py_modules = get_py_sources()
else:
    py_modules = []


ext_modules = []
for name in ["pyhelm_eos", "ic", "create_ics", "arepo_vis", "pyeos", "pyopal_eos", "pysph"]:
    ext_modules.append(
        Extension(name,
                  include_dirs=include_dirs,
                  libraries=libraries,
                  library_dirs=library_dirs,
                  define_macros=define_macros,
                  extra_compile_args=['-fopenmp'],
                  extra_link_args=['-fopenmp'],
                  sources=sources)
    )

setup(name="AREPO Helper Library",
      version="1.0",
      description="Various functions and class definitions to help create and analyse AREPO simulations.",
      author="Uri Pierre Burmester",
      author_email="uri.burmester@anu.edu.au",
      py_modules=py_modules,
      ext_modules=ext_modules)
