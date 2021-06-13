from numpy.distutils.core import setup, Extension
import numpy as np
import os

incl_dirs = ['/usr/include/python3.8/', '/usr/include/hdf5/serial/',
             os.path.join(np.get_include(), 'numpy')]
libs_dirs = ['/usr/lib/python3.8/']
libs = ['gsl', 'gslcblas', 'm', 'armadillo']
define_macros = [('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')]
c_libs_dir = 'arepo_helper/libs/'
py_dir = 'arepo_helper/'

pyhelm_eos = Extension('pyhelm_eos',
                       include_dirs=incl_dirs,
                       libraries=libs,
                       library_dirs=libs_dirs,
                       define_macros=define_macros,
                       sources=[c_libs_dir + 'pyhelm_eos.cpp',
                                c_libs_dir + 'helm_eos.cpp'])

ic = Extension('ic',
               include_dirs=incl_dirs,
               libraries=libs,
               library_dirs=libs_dirs,
               define_macros=define_macros,
               sources=[c_libs_dir + 'ic.cpp',
                        c_libs_dir + 'pyhelm_eos.cpp',
                        c_libs_dir + 'helm_eos.cpp'])

create_ics = Extension('create_ics',
                       include_dirs=incl_dirs,
                       libraries=libs,
                       library_dirs=libs_dirs,
                       define_macros=define_macros,
                       sources=[c_libs_dir + 'create_ics.cpp',
                                c_libs_dir + 'pyhelm_eos.cpp',
                                c_libs_dir + 'helm_eos.cpp'])

arepo_pcolor = Extension('arepo_pcolor',
                         include_dirs=incl_dirs,
                         libraries=libs,
                         library_dirs=libs_dirs,
                         define_macros=define_macros,
                         sources=[c_libs_dir + 'make_pcolor.cpp',
                                  c_libs_dir + 'sph.cpp',
                                  c_libs_dir + 'ic.cpp',
                                  c_libs_dir + 'pyhelm_eos.cpp',
                                  c_libs_dir + 'helm_eos.cpp'])

arepo_radial = Extension('arepo_radial',
                         include_dirs=incl_dirs,
                         libraries=libs,
                         library_dirs=libs_dirs,
                         define_macros=define_macros,
                         sources=[c_libs_dir + 'make_radial.cpp',
                                  c_libs_dir + 'sph.cpp',
                                  c_libs_dir + 'ic.cpp',
                                  c_libs_dir + 'pyhelm_eos.cpp',
                                  c_libs_dir + 'helm_eos.cpp'])

setup(name='AREPO Helper Library',
      version='1,0',
      description='Various functions and class definitions to ',
      author='Uri Pierre Burmester',
      author_email='uri.burmester@anu.edu.au',
      py_modules=[
          py_dir + 'scatter_2d',
          py_dir + 'abstract_plot',
          py_dir + 'analysis',
          py_dir + 'group_animation',
          py_dir + 'h5_file',
          py_dir + 'arepo_helper',
          py_dir + 'ics',
          py_dir + 'names',
          py_dir + 'pcolor_plot',
          py_dir + 'plot_manager',
          py_dir + 'plot_options',
          py_dir + 'radial_plot',
          py_dir + 'run',
          py_dir + 'snapshot',
          py_dir + 'species',
          py_dir + 'utilities',
          py_dir + 'const',
          py_dir + 'wd_utils',
          py_dir + 'wdec_results'
      ],
      ext_modules=[
          pyhelm_eos,
          ic,
          create_ics,
          arepo_pcolor,
          arepo_radial])
