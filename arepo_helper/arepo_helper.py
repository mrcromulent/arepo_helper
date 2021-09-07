from plot_manager import quick_radial, make_radial_time_series, quick_pcolor
from abstract_plot import RadialPlotOptions, RadialPlot
from utilities import SuppressStdout
from analysis import ArepoAnalyser
import matplotlib.pyplot as plt
from utilities import sci
from run import ArepoRun
from names import n
import numpy as np


def main():
    """Main script."""
    ar = ArepoRun.from_directory("/home/pierre/Desktop/output")

    quant_name = n.TEMPERATURE
    for i, s in enumerate(ar.snapshots):
        qp = quick_pcolor(s, quant_name)
        qp.save()
        qr = quick_radial(s, quant_name)
        qr.save()
    make_radial_time_series(ar, n.DENSITY)

    # from definitions import DATA_DIR
    # from astro_utils import rho_min, rho_max, e_min, e_max, alpha, beta
    # from species import ArepoSpeciesList
    # from matplotlib.colors import SymLogNorm
    # from matplotlib.cm import ScalarMappable
    # import pickle as pk
    # from names import n
    # import matplotlib.pyplot as plt
    # from pyhelm_eos import loadhelm_eos
    # from wdec_results import WDECResults
    # import numpy as np
    # import ic
    # import os
    # from ic import create_wdHe, create_wd
    #
    # helm_file = os.path.join("/home/pierre/Downloads/misc/helmholtz/helm_table.dat")
    # species_file = os.path.join(DATA_DIR, "eostable/species05.txt")
    # rho = 8e5
    # temp = 5e5
    # xnuc = [1.0, 0.0, 0.0, 0.0, 0.0]
    # eos = loadhelm_eos(helm_file, species_file, True)
    # test = create_wdHe(eos, rho, temp)
    # test2 = create_wd(eos, rho, temp, xnuc)
    # # print(test)
    #
    # fig, ax = plt.subplots()
    # # ax.plot(test['r'], test['p'])
    # # ax.plot(test2['Radius'], test2[n.PRESSURE])
    #
    # ax.plot(test['r'], test['u'])
    # ax.plot(test2['Radius'], test2[n.INTERNALENERGY])
    #
    # # ax.plot(test['r'], test['rho'])
    # # ax.plot(test2['Radius'], test2[n.DENSITY])
    #
    # ax.set_yscale('log')
    # fig.show()


if __name__ == "__main__":
    main()
