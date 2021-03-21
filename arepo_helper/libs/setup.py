from distutils.core import setup, Extension
import numpy as np
import os


# the first path can probably be removed because of the automatic detection below
incl_dirs = [
    '/home/pierre/PycharmProjects/arepo_helper/lib64/python3.8/site-packages/numpy/core/include/numpy',
    'home/pierre/PycharmProjects/arepo_helper/include']
libs_dirs = [
    'home/pierre/PycharmProjects/arepo_helper/lib']

# find numpy include path automatically
incl_dirs       += [os.path.join(np.get_include(), 'numpy')]
define_macros   = [("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")]
libs            = ['gsl', 'gslcblas', 'm']

grid = Extension('grid',
                 include_dirs=incl_dirs,
                 libraries=libs,
                 library_dirs=libs_dirs,
                 define_macros=define_macros,
                 sources=['grid.c'])

pyeos = Extension('pyeos',
                  include_dirs=incl_dirs,
                  libraries=libs,
                  library_dirs=libs_dirs,
                  define_macros=define_macros,
                  sources=['pyeos.c', 'eos.c'])

pyhelm_eos = Extension('pyhelm_eos',
                       include_dirs=incl_dirs,
                       libraries=libs,
                       library_dirs=libs_dirs,
                       define_macros=define_macros,
                       sources=['pyhelm_eos.c', 'helm_eos.c'])

pyopal_eos = Extension('pyopal_eos',
                       include_dirs=incl_dirs,
                       libraries=libs,
                       library_dirs=libs_dirs,
                       define_macros=define_macros,
                       sources=['pyopal_eos.c', 'opal_eos.c'])

pysph = Extension('pysph',
                  include_dirs=incl_dirs,
                  libraries=libs,
                  library_dirs=libs_dirs,
                  define_macros=define_macros,
                  sources=['pysph.c', 'sph.c'])

calcGrid = Extension('calcGrid',
                     include_dirs=incl_dirs,
                     libraries=libs,
                     library_dirs=libs_dirs,
                     define_macros=define_macros,
                     sources=['calcGrid.c', 'sph.c'])

ic = Extension('ic',
               include_dirs=incl_dirs,
               libraries=libs,
               library_dirs=libs_dirs,
               define_macros=define_macros,
               sources=['ic.c', 'pyhelm_eos.c', 'helm_eos.c'])

rgadget = Extension('rgadget',
                    include_dirs=incl_dirs,
                    libraries=libs,
                    library_dirs=libs_dirs,
                    define_macros=define_macros,
                    sources=['rgadget.c', 'gadgetSnap.c'])

boundmass = Extension('boundmass',
                      include_dirs=incl_dirs,
                      libraries=libs,
                      library_dirs=libs_dirs,
                      define_macros=define_macros,
                      sources=['boundmass.c', 'gadgetSnap.c'])

createICs = Extension('createICs',
                      include_dirs=incl_dirs,
                      libraries=libs,
                      library_dirs=libs_dirs,
                      define_macros=define_macros,
                      sources=['createICs.c', 'gadgetSnap.c', 'pyhelm_eos.c', 'helm_eos.c', 'rgadget.c'])

setup(name='Python script lib',
      version='1.0',
      description='Scripts to work with GADGET/LEAFS input & output',
      author='Ruediger',
      author_email='ruediger.pakmor@h-its.org',
      ext_modules=[grid, pyeos, pyhelm_eos, pyopal_eos, pysph, calcGrid, ic, rgadget, boundmass, createICs])
