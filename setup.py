from numpy.distutils.core import setup, Extension
import numpy as np
import os

incl_dirs   = ['/usr/include/python3.8/']
libs_dirs   = ['/usr/lib/python3.8/']
incl_dirs   += [os.path.join(np.get_include(), 'numpy')]
libs        = ['gsl', 'gslcblas', 'm']
defines     = [('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')]

grid = Extension('grid',
                 define_macros=defines,
                 include_dirs=incl_dirs,
                 libraries=libs,
                 library_dirs=libs_dirs,
                 sources=['libs/grid.c'])

pyeos = Extension('pyeos',
                  define_macros=defines,
                  include_dirs=incl_dirs,
                  libraries=libs,
                  library_dirs=libs_dirs,
                  sources=['libs/pyeos.c', 'libs/eos.c'])

pyhelm_eos = Extension('pyhelm_eos',
                       define_macros=defines,
                       include_dirs=incl_dirs,
                       libraries=libs,
                       library_dirs=libs_dirs,
                       sources=['libs/pyhelm_eos.c', 'libs/helm_eos.c'])

pyopal_eos = Extension('pyopal_eos',
                       define_macros=defines,
                       include_dirs=incl_dirs,
                       libraries=libs,
                       library_dirs=libs_dirs,
                       sources=['libs/pyopal_eos.c', 'libs/opal_eos.c'])

pysph = Extension('pysph',
                  define_macros=defines,
                  include_dirs=incl_dirs,
                  libraries=libs,
                  library_dirs=libs_dirs,
                  extra_compile_args=['-fopenmp'],
                  extra_link_args=['-fopenmp'],
                  sources=['libs/pysph.c', 'libs/sph.c'])

calcGrid = Extension('calcGrid',
                     define_macros=defines,
                     include_dirs=incl_dirs,
                     libraries=libs,
                     library_dirs=libs_dirs,
                     extra_compile_args=['-fopenmp'],
                     extra_link_args=['-fopenmp'],
                     sources=['libs/calcGrid.c', 'libs/sph.c'])

ic = Extension('ic',
               define_macros=defines,
               include_dirs=incl_dirs,
               libraries=libs,
               library_dirs=libs_dirs,
               sources=['libs/ic.c', 'libs/pyhelm_eos.c', 'libs/helm_eos.c'])

rgadget = Extension('rgadget',
                    define_macros=defines,
                    include_dirs=incl_dirs,
                    libraries=libs,
                    library_dirs=libs_dirs,
                    sources=['libs/rgadget.c', 'libs/gadgetSnap.c'])

boundmass = Extension('boundmass',
                      define_macros=defines,
                      include_dirs=incl_dirs,
                      libraries=libs,
                      library_dirs=libs_dirs,
                      sources=['libs/boundmass.c', 'libs/gadgetSnap.c'])

createICs = Extension('createICs',
                      define_macros=defines,
                      include_dirs=incl_dirs,
                      libraries=libs,
                      library_dirs=libs_dirs,
                      sources=['libs/createICs.c', 'libs/gadgetSnap.c', 'libs/pyhelm_eos.c', 'libs/helm_eos.c',
                               'libs/rgadget.c'])

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
      py_modules=['arepo_2d_scatter',
                  'arepo_abstract_plot',
                  'arepo_analysis',
                  'arepo_group_animation',
                  'arepo_group_frame',
                  'arepo_h5_file',
                  'arepo_helper',
                  'arepo_ics',
                  'arepo_names',
                  'arepo_pcolor_plot',
                  'arepo_plot_manager',
                  'arepo_plot_options',
                  'arepo_radial_plot',
                  'arepo_run',
                  'arepo_snapshot',
                  'arepo_species',
                  'arepo_utilities',
                  'const',
                  'wd_utils'],
      install_requires=['numpy==1.17.5'],
      ext_modules=[grid, pyeos, pyhelm_eos, pyopal_eos, pysph, calcGrid, ic, rgadget, boundmass, createICs],
      # opalopacities],
      data_files=[('eostable', ['eostable/EOS5_data', 'eostable/GN93hz',
                                'eostable/helm_table.dat', 'eostable/species_star.txt',
                                'eostable/species05.txt', 'eostable/species_wd_co.txt',
                                'eostable/decaydata', 'eostable/species384.txt',
                                'eostable/data_logR', 'eostable/data_logT',
                                'eostable/data_logkappa', 'eostable/leveldata'])]
      )
