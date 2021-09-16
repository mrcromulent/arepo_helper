from abstract_plot import AbstractPlot, get_plotter_func_from_plot_options
import matplotlib.animation as ani
import matplotlib.pyplot as plt
from datetime import datetime
from tqdm import tqdm
import numpy as np
import utilities
import os


class GroupAnimation:
    """Class to create a group of animations in a grid"""
    plots: np.ndarray = None
    pbar: tqdm = None
    animation: ani.FuncAnimation = None

    def __init__(self,
                 plot_options_array: np.ndarray,
                 fps: int = 5) -> None:
        """Constructor

        :param plot_options_array: Array of PlotOptions objects, where axis 2 indicates the passive of time
        :type plot_options_array: np.ndarray
        :param fps: The frames-per-second of the animation
        :type fps: int

        :return: Group animation object
        :rtype: GroupAnimation
        """

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

    def animate(self) -> None:
        """Creates the animation"""
        self.pbar = tqdm(total=self.num_frames)

        self.init()
        self.fig.tight_layout()
        self.animation = ani.FuncAnimation(self.fig,
                                           self.animate_frame,
                                           frames=self.num_frames,
                                           repeat=True,
                                           init_func=utilities.dummy)

    def init(self) -> None:
        """Initialises all plots to be animated"""

        for i in range(self.nrows):
            for j in range(self.ncols):
                po = self.poa[i, j, 0]
                plotter = get_plotter_func_from_plot_options(po)
                self.plots[i, j] = plotter(po, figure=self.fig, ax=self.ax[i, j])

    def animate_frame(self,
                      t: int) -> np.ndarray:
        """Updates the progress bar and animates a single frame

        :param t: The frame number to be animated
        :type t: int

        :return: Array of AbstractPlots
        :rtype: np.ndarray
        """

        # Update progress bar
        self.pbar.update(1)

        # Update all plots using their update methods
        for i in range(self.nrows):
            for j in range(self.ncols):
                po = self.poa[i, j, t]
                self.plots[i, j].update_plot_from_plot_options(po)

        return self.plots

    def save(self,
             filename: str = None) -> None:
        """Saves the animation to the file specified by filename

        :param filename: Location to save file to
        :type filename: str
        """

        if filename is None:
            dt = datetime.now().strftime("%Y%m%d-%H%M%S")
            fn_stem = f"{type(self).__name__}-{dt}.mp4"
            dir_name = os.path.dirname(self.poa[0, 0, 0].ar.snapshots[0].filename)
            filename = os.path.join(dir_name, fn_stem)

        self.animation.save(filename, fps=self.fps, metadata={"Comment": f"Animation"})
