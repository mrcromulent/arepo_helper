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
    ar = ArepoRun.from_directory("/home/pierre/Desktop/output_new_relaxobject_c")

    quant_name = n.TEMPERATURE
    for i, s in enumerate(ar.snapshots):
        qp = quick_pcolor(s, quant_name)
        qp.save()


if __name__ == "__main__":
    main()
