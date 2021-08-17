from abstract_plot import AbstractPlot, get_plotter_func_from_plot_options
import matplotlib.animation as ani
import matplotlib.pyplot as plt
from datetime import datetime
from tqdm import tqdm
import numpy as np
import utilities
import os


class GroupAnimation(object):
    plots = []
    pbar = None
    animation = None

    def __init__(self, plot_options_array, fps=5):

        self.poa    = plot_options_array
        self.fps    = fps
        self.nrows = np.size(plot_options_array, axis=0)
        self.ncols = np.size(plot_options_array, axis=1)
        self.num_frames = np.size(plot_options_array, axis=2)
        self.plots = np.empty((self.nrows, self.ncols), dtype=AbstractPlot)

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
                                           init_func=utilities.dummy)

    def init(self):

        for i in range(self.nrows):
            for j in range(self.ncols):
                po = self.poa[i, j, 0]
                plotter = get_plotter_func_from_plot_options(po)
                self.plots[i, j] = plotter(po, figure=self.fig, ax=self.ax[i, j])

    def animate_frame(self, t):

        # Update progress bar
        self.pbar.update(1)

        # Update all plots using their update methods
        for i in range(self.nrows):
            for j in range(self.ncols):
                po = self.poa[i, j, t]
                self.plots[i, j].update_plot_from_plot_options(po)

        return self.plots

    def save(self, filename=None):

        if filename is None:
            dt = datetime.now().strftime("%Y%m%d-%H%M%S")
            fn_stem = f"{type(self).__name__}-{dt}.mp4"
            dir_name = os.path.dirname(self.poa[0, 0, 0].ar.snapshots[0].filename)
            filename = os.path.join(dir_name, fn_stem)

        self.animation.save(filename, fps=self.fps, metadata={"Comment": f"Animation"})
