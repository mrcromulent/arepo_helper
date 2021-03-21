from numpy.distutils.core import setup, Extension
import numpy as np
import os

incl_dirs   = ['/usr/include/python3.8/']
libs_dirs   = ['/usr/lib/python3.8/']
incl_dirs   += [os.path.join(np.get_include(), 'numpy')]
libs        = ['gsl', 'gslcblas', 'm']
defines     = [('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')]
data_dir    = 'arepo_helper/data/eostable/'
arepo_libs  = 'arepo_helper/libs/'

grid = Extension('grid',
                 define_macros=defines,
                 include_dirs=incl_dirs,
                 libraries=libs,
                 library_dirs=libs_dirs,
                 sources=[arepo_libs+'grid.c'])

pyeos = Extension('pyeos',
                  define_macros=defines,
                  include_dirs=incl_dirs,
                  libraries=libs,
                  library_dirs=libs_dirs,
                  sources=[arepo_libs+'pyeos.c', arepo_libs+'eos.c'])

pyhelm_eos = Extension('pyhelm_eos',
                       define_macros=defines,
                       include_dirs=incl_dirs,
                       libraries=libs,
                       library_dirs=libs_dirs,
                       sources=[arepo_libs+'pyhelm_eos.c', arepo_libs+'helm_eos.c'])

pyopal_eos = Extension('pyopal_eos',
                       define_macros=defines,
                       include_dirs=incl_dirs,
                       libraries=libs,
                       library_dirs=libs_dirs,
                       sources=[arepo_libs+'pyopal_eos.c', arepo_libs+'opal_eos.c'])

pysph = Extension('pysph',
                  define_macros=defines,
                  include_dirs=incl_dirs,
                  libraries=libs,
                  library_dirs=libs_dirs,
                  extra_compile_args=['-fopenmp'],
                  extra_link_args=['-fopenmp'],
                  sources=[arepo_libs+'pysph.c', arepo_libs+'sph.c'])

calcGrid = Extension('calcGrid',
                     define_macros=defines,
                     include_dirs=incl_dirs,
                     libraries=libs,
                     library_dirs=libs_dirs,
                     extra_compile_args=['-fopenmp'],
                     extra_link_args=['-fopenmp'],
                     sources=[arepo_libs+'calcGrid.c', arepo_libs+'sph.c'])

ic = Extension('ic',
               define_macros=defines,
               include_dirs=incl_dirs,
               libraries=libs,
               library_dirs=libs_dirs,
               sources=[arepo_libs+'ic.c', arepo_libs+'pyhelm_eos.c', arepo_libs+'helm_eos.c'])

rgadget = Extension('rgadget',
                    define_macros=defines,
                    include_dirs=incl_dirs,
                    libraries=libs,
                    library_dirs=libs_dirs,
                    sources=[arepo_libs+'rgadget.c', arepo_libs+'gadgetSnap.c'])

boundmass = Extension('boundmass',
                      define_macros=defines,
                      include_dirs=incl_dirs,
                      libraries=libs,
                      library_dirs=libs_dirs,
                      sources=[arepo_libs+'boundmass.c', arepo_libs+'gadgetSnap.c'])

createICs = Extension('createICs',
                      define_macros=defines,
                      include_dirs=incl_dirs,
                      libraries=libs,
                      library_dirs=libs_dirs,
                      sources=[arepo_libs+'createICs.c',
                               arepo_libs+'gadgetSnap.c',
                               arepo_libs+'pyhelm_eos.c',
                               arepo_libs+'helm_eos.c',
                               arepo_libs+'rgadget.c'])

# opalopacities = Extension('opalopacities',
#                           define_macros=defines,
#                           include_dirs=incl_dirs,
#                           libraries=libs,
#                           library_dirs=libs_dirs,
#                           sources=['libs/opacities/opacities.pyf', 'libs/opacities/xztrin21.f'])

setup(name='AREPO Helper Library',
      version='0.1',
      description='Various functions and class definitions to ',
      author='Uri Pierre Burmester',
      author_email='uri.burmester@anu.edu.au',
      py_modules=['scatter_2d',
                  'abstract_plot',
                  'analysis',
                  'group_animation',
                  'group_frame',
                  'h5_file',
                  'arepo_helper',
                  'ics',
                  'names',
                  'pcolor_plot',
                  'plot_manager',
                  'plot_options',
                  'radial_plot',
                  'run',
                  'snapshot',
                  'species',
                  'utilities',
                  'const',
                  'wd_utils',
                  'wdec_results'],
      install_requires=['numpy==1.17.5'],
      ext_modules=[grid, pyeos, pyhelm_eos, pyopal_eos, pysph, calcGrid, ic, rgadget, boundmass, createICs],
      # opalopacities],
      data_files=[('eostable', [data_dir+'EOS5_data', data_dir+'GN93hz',
                                data_dir+'helm_table.dat', data_dir+'species_star.txt',
                                data_dir+'species05.txt', data_dir+'species_wd_co.txt',
                                data_dir+'decaydata', data_dir+'species384.txt',
                                data_dir+'data_logR', data_dir+'data_logT',
                                data_dir+'data_logkappa', data_dir+'leveldata'])]
      )
