from setuptools import setup, Extension
from sysconfig import get_paths
import numpy as np
import subprocess
import glob
import sys
import os


def get_cpp_sources():
    return [file for file in sorted(glob.glob("src/arepo_helper/cpp_libs/src/*.cpp"))]


def get_cpp_headers():
    return [file for file in sorted(glob.glob("src/arepo_helper/cpp_libs/src/*.h"))]


HDF5_dir = os.environ.get('HDF5_DIR')
HDF5_incdir = os.environ.get('HDF5_INCDIR')
HDF5_libdir = os.environ.get('HDF5_LIBDIR')

try:
    HAS_PKG_CONFIG = subprocess.call(['pkg-config', '--libs', 'hdf5'],
                                     stdout=subprocess.PIPE) == 0
except OSError:
    HAS_PKG_CONFIG = False

open_kwargs = {'encoding': 'utf-8'}


def check_hdf5version(hdf5_includedir):
    try:
        f = open(os.path.join(hdf5_includedir, 'H5public.h'), **open_kwargs)
    except IOError:
        return None
    hdf5_version = None
    for line in f:
        if line.startswith('#define H5_VERS_INFO'):
            hdf5_version = line.split('"')[1]
    return hdf5_version


def get_hdf5_version(direc):
    # check to see if hdf5 headers in direc, return version number or None
    sys.stdout.write('checking %s ...\n' % direc)
    hdf5_version = check_hdf5version(direc)
    if hdf5_version is None:
        sys.stdout.write('hdf5 headers not found in %s\n' % direc)
        return None
    else:
        sys.stdout.write('%s headers found in %s\n' % (hdf5_version, direc))
        return hdf5_version


def _populate_hdf5_info(dirstosearch, inc_dirs, libs, lib_dirs):
    global HDF5_incdir, HDF5_dir, HDF5_libdir

    nohdf5dirs = HDF5_incdir is None and HDF5_libdir is None and HDF5_dir is None
    if HAS_PKG_CONFIG and nohdf5dirs:
        # if HDF5 dirs not specified, and pkg-config available, use it
        dep = subprocess.Popen(['pkg-config', '--cflags', 'hdf5'],
                               stdout=subprocess.PIPE).communicate()[0]
        inc_dirs.extend([str(i[2:].decode()) for i in dep.split() if
                         i[0:2].decode() == '-I'])
        dep = subprocess.Popen(['pkg-config', '--libs', 'hdf5'],
                               stdout=subprocess.PIPE).communicate()[0]
        libs.extend(
            [str(l[2:].decode()) for l in dep.split() if l[0:2].decode() == '-l'])
        lib_dirs.extend(
            [str(l[2:].decode()) for l in dep.split() if l[0:2].decode() == '-L'])
        dep = subprocess.Popen(['pkg-config', '--cflags', 'hdf5'],
                               stdout=subprocess.PIPE).communicate()[0]
        inc_dirs.extend(
            [str(i[2:].decode()) for i in dep.split() if i[0:2].decode() == '-I'])
    else:
        if HDF5_incdir is None and HDF5_dir is None:
            sys.stdout.write("""
    HDF5_DIR environment variable not set, checking some standard locations ..\n""")
            for direc in dirstosearch:
                hdf5_version = get_hdf5_version(os.path.join(direc, 'include'))
                if hdf5_version is None:
                    continue
                else:
                    HDF5_dir = direc
                    HDF5_incdir = os.path.join(direc, 'include')
                    sys.stdout.write('%s found in %s\n' % (hdf5_version, HDF5_dir))
                    break
            if HDF5_dir is None:
                raise ValueError('did not find HDF5 headers')
        else:
            if HDF5_incdir is None:
                HDF5_incdir = os.path.join(HDF5_dir, 'include')
            hdf5_version = get_hdf5_version(HDF5_incdir)
            if hdf5_version is None:
                raise ValueError('did not find HDF5 headers in %s' % HDF5_incdir)
            else:
                sys.stdout.write('%s found in %s\n' % (hdf5_version, HDF5_dir))

        if HDF5_libdir is None and HDF5_dir is not None:
            HDF5_libdir = os.path.join(HDF5_dir, 'lib')

        if HDF5_libdir is not None:
            lib_dirs.append(HDF5_libdir)
        if HDF5_incdir is not None:
            inc_dirs.append(HDF5_incdir)

        libs.extend(['hdf5_hl', 'hdf5'])

    def _populate_armadillo_info():
        pass


dirstosearch = [os.path.expanduser('~'), '/usr/local', '/sw', '/opt', '/opt/local', '/usr']

lib_dirs = []
inc_dirs = []
libs = []

# _populate_hdf5_info will use HDF5_dir, HDF5_libdir and HDF5_incdir if they are set.
# otherwise pkg-config will be tried, and if that fails, dirstosearch will be searched.
_populate_hdf5_info(dirstosearch, inc_dirs, libs, lib_dirs)

# Getting Python and numpy src dirs
info = get_paths()
python_incl = info["include"]
python_lib = info["stdlib"]
numpy_incl = os.path.join(np.get_include(), "numpy")

# Source files and definitions for this thing
include_dirs = [python_incl, numpy_incl, "/home/149/ub0692/include/", inc_dirs[0], "src/arepo_helper/cpp_libs/include/"]
library_dirs = [python_lib, "/home/149/ub0692/lib64", lib_dirs[0]]
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
