from setuptools import setup, Extension
from sysconfig import get_paths
import numpy as np
import glob
import os


def get_cpp_sources():
    return [os.path.abspath(file) for file in sorted(glob.glob("./src/cpp_libs/src/*.cpp"))]


def get_cpp_headers():
    return [os.path.abspath(file) for file in sorted(glob.glob("./src/cpp_libs/src/*.h"))]

# NCI
# install_py_modules = True
# armadillo_incl = "/home/149/ub0692/include/"
# armadillo_lib = "/home/149/ub0692/lib64"
# hdf5_incl = "/apps/hdf5/1.10.5/include"
# hdf5_lib_dir = "/apps/hdf5/1.10.5/lib"

# Home
install_py_modules = True
hdf5_incl = "/usr/include/hdf5/serial/"
hdf5_lib_dir = "/usr/lib/x86_64-linux-gnu/hdf5/serial/"
armadillo_incl = "/"
armadillo_lib = "/"

# Getting Python and numpy src dirs
info = get_paths()
python_incl = info["include"]
python_lib = info["stdlib"]
numpy_incl = os.path.join(np.get_include(), "numpy")

# Source files and definitions for this thing
include_dirs = [python_incl, numpy_incl, armadillo_incl, hdf5_incl, os.path.abspath("./src/cpp_libs/src/")]
library_dirs = [python_lib, armadillo_lib, hdf5_lib_dir]
libraries = ["gsl", "gslcblas", "m", "armadillo", "hdf5"]
define_macros = [("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION"), ('H5_BUILT_AS_DYNAMIC_LIB', True)]
sources = get_cpp_sources()

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

setup(ext_modules=ext_modules)
