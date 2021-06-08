from distutils.core import setup, Extension
import numpy as np
import os


# the first path can probably be removed because of the automatic detection below
incl_dirs = [
    '/home/pierre/PycharmProjects/arepo_helper/lib64/python3.8/site-packages/numpy/core/include/numpy',
    '/home/pierre/PycharmProjects/arepo_helper/include',
    '/home/pierre/PycharmProjects/arepo_helper/arepo_helper/libs',
    '/usr/include/hdf5/serial/']
libs_dirs = ['/home/pierre/PycharmProjects/arepo_helper/lib/']

# find numpy include path automatically
incl_dirs       += [os.path.join(np.get_include(), 'numpy')]
define_macros   = [("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")]
libs            = ['gsl', 'gslcblas', 'm', 'armadillo']

pyhelm_eos = Extension('pyhelm_eos',
                       include_dirs=incl_dirs,
                       libraries=libs,
                       library_dirs=libs_dirs,
                       define_macros=define_macros,
                       sources=['pyhelm_eos.cpp', 'helm_eos.cpp'])

calcGrid = Extension('calcGrid',
                     include_dirs=incl_dirs,
                     libraries=libs,
                     library_dirs=libs_dirs,
                     define_macros=define_macros,
                     sources=['calcGrid.cpp', 'sph.cpp'])

ic = Extension('ic',
               include_dirs=incl_dirs,
               libraries=libs,
               library_dirs=libs_dirs,
               define_macros=define_macros,
               sources=['ic.cpp', 'pyhelm_eos.cpp', 'helm_eos.cpp'])

createICs = Extension('createICs',
                      include_dirs=incl_dirs,
                      libraries=libs,
                      library_dirs=libs_dirs,
                      define_macros=define_macros,
                      sources=['createICs.cpp', 'pyhelm_eos.cpp', 'helm_eos.cpp'])

makePColorPierre = Extension('pcolor_pierre',
                       include_dirs=incl_dirs,
                       libraries=libs,
                       library_dirs=libs_dirs,
                       define_macros=define_macros,
                       sources=['make_pcolor.cpp', 'sph.cpp', 'ic.cpp', 'pyhelm_eos.cpp', 'helm_eos.cpp'])


setup(name='Python script lib',
      version='1.0',
      description='Scripts to work with GADGET/LEAFS input & output',
      author='Ruediger',
      author_email='ruediger.pakmor@h-its.org',
      ext_modules=[
          pyhelm_eos,
          calcGrid,
          ic,
          createICs,
          makePColorPierre])
