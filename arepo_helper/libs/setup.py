from distutils.core import setup, Extension
import numpy as np
import os

# the first path can probably be removed because of the automatic detection below
incl_dirs = [
    '/home/pierre/PycharmProjects/arepo_helper/lib64/python3.8/site-packages/numpy/core/include/numpy',
    # '/home/pierre/PycharmProjects/arepo_helper/include',
    '/home/pierre/PycharmProjects/arepo_helper/arepo_helper/libs',
    '/usr/include/hdf5/serial/']
libs_dirs = ['/home/pierre/PycharmProjects/arepo_helper/lib/']

# find numpy include path automatically
incl_dirs += [os.path.join(np.get_include(), 'numpy')]
define_macros = [("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")]
libs = ['gsl', 'gslcblas', 'm', 'armadillo']

pyhelm_eos = Extension('pyhelm_eos',
                       include_dirs=incl_dirs,
                       libraries=libs,
                       library_dirs=libs_dirs,
                       define_macros=define_macros,
                       sources=['pyhelm_eos.cpp', 'helm_eos.cpp'])

ic = Extension('ic',
               include_dirs=incl_dirs,
               libraries=libs,
               library_dirs=libs_dirs,
               define_macros=define_macros,
               sources=['ic.cpp', 'pyhelm_eos.cpp', 'helm_eos.cpp'])

create_ics = Extension('create_ics',
                       include_dirs=incl_dirs,
                       libraries=libs,
                       library_dirs=libs_dirs,
                       define_macros=define_macros,
                       sources=['create_ics.cpp', 'pyhelm_eos.cpp', 'helm_eos.cpp'])

pcolor_pierre = Extension('pcolor_pierre',
                          include_dirs=incl_dirs,
                          libraries=libs,
                          library_dirs=libs_dirs,
                          define_macros=define_macros,
                          sources=['make_pcolor.cpp', 'sph.cpp', 'ic.cpp', 'pyhelm_eos.cpp', 'helm_eos.cpp'])

radial_pierre = Extension('radial_pierre',
                          include_dirs=incl_dirs,
                          libraries=libs,
                          library_dirs=libs_dirs,
                          define_macros=define_macros,
                          sources=['make_radial.cpp', 'sph.cpp', 'ic.cpp', 'pyhelm_eos.cpp', 'helm_eos.cpp'])

setup(name='Arepo Helper Libs',
      version='1.0',
      description='Arepo Helper Libs',
      author='Uri Pierre Burmester',
      author_email='uri.burmester@anu.edu.au',
      ext_modules=[
          pyhelm_eos,
          ic,
          create_ics,
          pcolor_pierre,
          radial_pierre])
