from numpy.distutils.core import setup, Extension
import numpy as np
import os

incl_dirs   = ['/usr/include/python3.8/', '/usr/include/hdf5/serial/']
libs_dirs   = ['/usr/lib/python3.8/']
incl_dirs   += [os.path.join(np.get_include(), 'numpy')]
libs        = ['gsl', 'gslcblas', 'm', 'armadillo']
defines     = [('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')]
data_dir    = 'arepo_helper/data/eostable/'
arepo_libs  = 'arepo_helper/libs/'
py_dir      = 'arepo_helper/'

pyhelm_eos = Extension('pyhelm_eos',
                       include_dirs=incl_dirs,
                       libraries=libs,
                       library_dirs=libs_dirs,
                       define_macros=defines,
                       sources=[arepo_libs+'pyhelm_eos.cpp', arepo_libs+'helm_eos.cpp'])

calcGrid = Extension('calcGrid',
                     include_dirs=incl_dirs,
                     libraries=libs,
                     library_dirs=libs_dirs,
                     define_macros=defines,
                     sources=[arepo_libs+'calcGrid.cpp', arepo_libs+'sph.cpp'])

ic = Extension('ic',
               include_dirs=incl_dirs,
               libraries=libs,
               library_dirs=libs_dirs,
               define_macros=defines,
               sources=[arepo_libs+'ic.cpp', arepo_libs+'pyhelm_eos.cpp', arepo_libs+'helm_eos.cpp'])

createICs = Extension('createICs',
                      include_dirs=incl_dirs,
                      libraries=libs,
                      library_dirs=libs_dirs,
                      define_macros=defines,
                      sources=[arepo_libs+'createICs.cpp', arepo_libs+'pyhelm_eos.cpp', arepo_libs+'helm_eos.cpp'])

createICs2 = Extension('createICs2',
                      include_dirs=incl_dirs,
                      libraries=libs,
                      library_dirs=libs_dirs,
                      define_macros=defines,
                      sources=[arepo_libs+'createICs.cpp', arepo_libs+'pyhelm_eos.cpp', arepo_libs+'helm_eos.cpp'])

makePColorPierre = Extension('pcolor_pierre',
                       include_dirs=incl_dirs,
                       libraries=libs,
                       library_dirs=libs_dirs,
                       define_macros=defines,
                       sources=[arepo_libs+'make_pcolor.cpp', arepo_libs+'sph.cpp', arepo_libs+'ic.cpp',
                                arepo_libs+'pyhelm_eos.cpp', arepo_libs+'helm_eos.cpp'])



setup(name='AREPO Helper Library',
      version='0.1',
      description='Various functions and class definitions to ',
      author='Uri Pierre Burmester',
      author_email='uri.burmester@anu.edu.au',
      py_modules=[
            py_dir+'scatter_2d',
            py_dir+'abstract_plot',
            py_dir+'analysis',
            py_dir+'group_animation',
            # py_dir+'group_frame',
            py_dir+'h5_file',
            py_dir+'arepo_helper',
            py_dir+'ics',
            py_dir+'names',
            py_dir+'pcolor_plot',
            py_dir+'plot_manager',
            py_dir+'plot_options',
            py_dir+'radial_plot',
            py_dir+'run',
            py_dir+'snapshot',
            py_dir+'species',
            py_dir+'utilities',
            py_dir+'const',
            py_dir+'wd_utils',
            py_dir+'wdec_results'
      ],
      install_requires=['numpy==1.17.5'],
      ext_modules=[
          pyhelm_eos,
          calcGrid,
          ic,
          createICs,
          createICs2,
          makePColorPierre],
      # data_files=[('eostable', [data_dir+'EOS5_data', data_dir+'GN93hz',
      #                           data_dir+'helm_table.dat', data_dir+'species_star.txt',
      #                           data_dir+'species05.txt', data_dir+'species_wd_co.txt',
      #                           data_dir+'decaydata', data_dir+'species384.txt',
      #                           data_dir+'data_logR', data_dir+'data_logT',
      #                           data_dir+'data_logkappa', data_dir+'leveldata'])]
      )
