from plot_options import ScatterPlotOptions, RadialPlotOptions, PColorPlotOptions
from pcolor_plot import PColorPlot
from radial_plot import RadialPlot
from scatter_2d import Scatter2D
import matplotlib.pyplot as plt
from names import n
import numpy as np


class GroupFrame(object):

    def __init__(self, plot_options_array):

        self.tvals = []
        self.orientations = []
        self.quantities = []

        if plot_options_array.ndim == 1:
            plot_options_array = np.expand_dims(plot_options_array, axis=1)

        for idx, po in np.ndenumerate(plot_options_array):
            self.tvals.append(po.t)
            self.orientations.append(po.orientation)
            self.quantities.append(po.quantity)

        self.nrows = np.size(plot_options_array, axis=0)
        self.ncols = np.size(plot_options_array, axis=1)

        self.plot_options_array = plot_options_array
        self.fig, self.ax = plt.subplots(nrows=self.nrows,
                                         ncols=self.ncols,
                                         figsize=(8 * self.ncols, 8 * self.nrows),
                                         squeeze=False)

        self.populate_plots()

    def populate_plots(self):

        for idx, po in np.ndenumerate(self.plot_options_array):
            # print(idx)
            ax = self.ax[idx]

            if isinstance(po, ScatterPlotOptions):
                Scatter2D(po, figure=self.fig, ax=ax)
            elif isinstance(po, PColorPlotOptions):
                PColorPlot(po, figure=self.fig, ax=ax)
            elif isinstance(po, RadialPlotOptions):
                RadialPlot(po, figure=self.fig, ax=ax)
            else:
                raise ValueError("Unknown plot options type")

        plt.tight_layout()

    def save(self, filename=None):
        if filename is None:
            filename = f"./{type(self).__name__}_{self.tvals}_{self.quantities}_{self.orientations}.png"
        self.fig.savefig(filename, dpi=300)


if __name__ == "__main__":
    from run import ArepoRun
    from analysis import ArepoAnalyser
    from plot_manager import PlotManager

    ar = ArepoRun.from_directory("./snapshots")
    aa = ArepoAnalyser(analysis_options={"inner_boxsize": 5e9})
    apm = PlotManager(ar, analyser=aa)

    for t in [0, 10, 20, 30, 40, 50, 80, 100, 130]:

        poa = np.transpose(np.array((
            [apm.compute_plot_options(t, n.DENSITY, ["x", "y"], "Radial", explicit_options={"logscale": True})],
            [apm.compute_plot_options(t, n.DENSITY, ["x", "y"], "Scatter", explicit_options={"log_cmap": True})],
            [apm.compute_plot_options(t, n.DENSITY, ["x", "y"], "PColor", explicit_options={"log_cmap": True})]
        )))

        test = GroupFrame(poa)
        test.save()
