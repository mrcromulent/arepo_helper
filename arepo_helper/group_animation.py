from plot_options import PlotOptions, ScatterPlotOptions, RadialPlotOptions, PColorPlotOptions
from abstract_plot import AbstractPlot
from pcolor_plot import PColorPlot
from radial_plot import RadialPlot
from scatter_2d import Scatter2D
import matplotlib.animation as ani
from utilities import dummy
import matplotlib.pyplot as plt
from names import n
from tqdm import tqdm
import numpy as np


class GroupAnimation(object):
    plots = []
    pbar = None
    animation = None

    def __init__(self, trange, plot_options_array, arepo_plot_manager):

        # plot_options_array = plot_options_array.reshape((1, 1))
        if plot_options_array.ndim == 1:
            plot_options_array = np.expand_dims(plot_options_array, axis=1)

        self.trange = trange
        self.apm = arepo_plot_manager
        self.num_frames = trange[1] - trange[0]

        # Initialise empty arrays for plot options and plot objects
        self.poa = np.empty((*plot_options_array.shape, self.num_frames), dtype=PlotOptions)
        self.plots = np.empty(plot_options_array.shape, dtype=AbstractPlot)

        for idx, po in np.ndenumerate(plot_options_array):
            for t, tval in enumerate(range(trange[0], trange[1])):
                q = po.quantity
                o = po.orientation

                self.poa[idx[0], idx[1], t] = self.apm.compute_plot_options(tval, q, o, po.plot_type,
                                                                            explicit_options=po.explicit_options)

        self.nrows = np.size(plot_options_array, axis=0)
        self.ncols = np.size(plot_options_array, axis=1)
        self.fig, self.ax = plt.subplots(nrows=self.nrows,
                                         ncols=self.ncols,
                                         figsize=(8 * self.ncols, 8 * self.nrows),
                                         squeeze=False)

    def animate(self):
        self.pbar = tqdm(total=self.num_frames)

        self.init()
        self.fig.tight_layout()
        self.animation = ani.FuncAnimation(self.fig,
                                           self.animate_frame,
                                           frames=self.num_frames,
                                           repeat=True,
                                           init_func=dummy)

    def init(self):

        for i in range(self.nrows):
            for j in range(self.ncols):

                ax = self.ax[i, j]
                po = self.poa[i, j, 0]

                if isinstance(po, ScatterPlotOptions):
                    self.plots[i, j] = Scatter2D(po, figure=self.fig, ax=ax)
                elif isinstance(po, PColorPlotOptions):
                    self.plots[i, j] = PColorPlot(po, figure=self.fig, ax=ax)
                elif isinstance(po, RadialPlotOptions):
                    self.plots[i, j] = RadialPlot(po, figure=self.fig, ax=ax)
                else:
                    raise ValueError("Unknown plot options type")

                # ax.set_aspect('equal')

    def animate_frame(self, t):

        self.pbar.update(1)

        for i in range(self.nrows):
            for j in range(self.ncols):
                po = self.poa[i, j, t]
                self.plots[i, j].update_plot_from_plot_options(po)

                # plt.tight_layout()

        return self.plots

    def save(self, filename=None):

        if filename is None:
            filename = f"./{type(self).__name__}_{self.trange}.mp4"

        metadata = {"Comment": f"Animation with : \n {str(self.apm)}"}
        # metadata = dict()
        self.animation.save(filename, fps=5, metadata=metadata)


if __name__ == "__main__":
    from run import ArepoRun
    from analysis import ArepoAnalyser
    from plot_manager import PlotManager

    ar = ArepoRun.from_directory("./snapshots")
    aa = ArepoAnalyser(analysis_options={"inner_boxsize": 5e9})
    apm = PlotManager(ar, analyser=aa)

    poa = np.transpose(np.array((
        [apm.compute_plot_options(0, n.PRESSURE, ["x", "y"], "Radial", explicit_options={"logscale": True})],
        [apm.compute_plot_options(0, n.PRESSURE, ["x", "y"], "Scatter", explicit_options={"log_cmap": True})],
        [apm.compute_plot_options(0, n.PRESSURE, ["x", "y"], "PColor", explicit_options={"log_cmap": True})]
    )))

    test = GroupAnimation([0, 130], poa, apm)
    test.animate()
    test.save()
