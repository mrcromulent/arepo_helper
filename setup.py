from numpy.distutils.core import setup, Extension
from definitions import ROOT_DIR
from sysconfig import get_paths
import numpy as np
import os

# Getting Python and numpy src dirs
info = get_paths()
python_incl = info["include"]
python_lib = info["stdlib"]
numpy_incl = os.path.join(np.get_include(), 'numpy')
c_src_dir = os.path.join(ROOT_DIR, 'arepo_helper/libs')
py_src_dir = os.path.join(ROOT_DIR, 'arepo_helper')

# NCI
# install_py_modules = True
# armadillo_incl = '/home/149/ub0692/include/'
# armadillo_lib = '/home/149/ub0692/lib64'
# hdf5_incl = '/apps/hdf5/1.10.5/include'

# Home
install_py_modules = False
hdf5_incl = '/usr/include/hdf5/serial/'
armadillo_incl = '/'
armadillo_lib = '/'

# Source files and definitions for this thing
incl_dirs = [python_incl, numpy_incl, armadillo_incl, hdf5_incl, c_src_dir]
libs_dirs = [python_lib, armadillo_lib]
libs = ['gsl', 'gslcblas', 'm', 'armadillo']
define_macros = [('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')]


if install_py_modules:
    py_modules = [
        os.path.join(py_src_dir, 'scatter_2d'),
        os.path.join(py_src_dir, 'abstract_plot'),
        os.path.join(py_src_dir, 'analysis'),
        os.path.join(py_src_dir, 'group_animation'),
        os.path.join(py_src_dir, 'h5_file'),
        os.path.join(py_src_dir, 'arepo_helper'),
        os.path.join(py_src_dir, 'ics'),
        os.path.join(py_src_dir, 'names'),
        os.path.join(py_src_dir, 'pcolor_plot'),
        os.path.join(py_src_dir, 'plot_manager'),
        os.path.join(py_src_dir, 'plot_options'),
        os.path.join(py_src_dir, 'radial_plot'),
        os.path.join(py_src_dir, 'run'),
        os.path.join(py_src_dir, 'snapshot'),
        os.path.join(py_src_dir, 'species'),
        os.path.join(py_src_dir, 'utilities'),
        os.path.join(py_src_dir, 'const'),
        os.path.join(py_src_dir, 'wd_utils'),
        os.path.join(py_src_dir, 'wdec_results')
    ]
else:
    py_modules = []

ext_modules = [
    Extension('pyhelm_eos',
              include_dirs=incl_dirs,
              libraries=libs,
              library_dirs=libs_dirs,
              define_macros=define_macros,
              sources=[
                  os.path.join(c_src_dir, 'pyhelm_eos.cpp'),
                  os.path.join(c_src_dir, 'helm_eos.cpp'),
                  os.path.join(c_src_dir, 'const.cpp')
              ]),
    Extension('ic',
              include_dirs=incl_dirs,
              libraries=libs,
              library_dirs=libs_dirs,
              define_macros=define_macros,
              sources=[
                  os.path.join(c_src_dir, 'ic.cpp'),
                  os.path.join(c_src_dir, 'pyhelm_eos.cpp'),
                  os.path.join(c_src_dir, 'helm_eos.cpp'),
                  os.path.join(c_src_dir, 'const.cpp')
              ]),
    Extension('create_ics',
              include_dirs=incl_dirs,
              libraries=libs,
              library_dirs=libs_dirs,
              define_macros=define_macros,
              sources=[
                  os.path.join(c_src_dir, 'create_ics.cpp'),
                  os.path.join(c_src_dir, 'pyhelm_eos.cpp'),
                  os.path.join(c_src_dir, 'helm_eos.cpp'),
                  os.path.join(c_src_dir, 'const.cpp')
              ]),
    Extension('arepo_vis',
              include_dirs=incl_dirs,
              libraries=libs,
              library_dirs=libs_dirs,
              define_macros=define_macros,
              sources=[
                  os.path.join(c_src_dir, 'visualise.cpp'),
                  os.path.join(c_src_dir, 'sph.cpp'),
                  os.path.join(c_src_dir, 'ic.cpp'),
                  os.path.join(c_src_dir, 'pyhelm_eos.cpp'),
                  os.path.join(c_src_dir, 'helm_eos.cpp'),
                  os.path.join(c_src_dir, 'const.cpp')
              ])
]

setup(name='AREPO Helper Library',
      version='1.0',
      description='Various functions and class definitions to help create and analyse AREPO simulations.',
      author='Uri Pierre Burmester',
      author_email='uri.burmester@anu.edu.au',
      py_modules=py_modules,
      ext_modules=ext_modules)
