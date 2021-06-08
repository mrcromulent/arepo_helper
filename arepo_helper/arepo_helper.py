"""
    # analysis_options = {"inner_boxsize": 1e10,
    #                     "cutoff_table": {n.TEMPERATURE: [0, 10e10]},
    #                     "highlights": {n.DENSITY: [2e6, 10e6]},
                          "slice_width": 0.002}

"""

from run import ArepoRun
from analysis import ArepoAnalyser
from plot_manager import PlotManager
# from group_frame import GroupFrame
from group_animation import GroupAnimation
from plot_options import PColorPlotOptions, RadialPlotOptions, ScatterPlotOptions
from scatter_2d import Scatter2D
from radial_plot import RadialPlot
from pcolor_plot import PColorPlot


class ArepoHelper(object):

    def __init__(self, arepo_plot_manager):
        self.apm = arepo_plot_manager

    def export_group_animation(self, ts, qs, os, plot_type, filename=None):

        # TODO: Fix array computation

        poa = self.apm.compute_plot_options_array([ts[0]], qs, os, plot_type)
        group_animation = GroupAnimation(ts, poa, self.apm)
        group_animation.animate()
        group_animation.save(filename)

    def export_animation(self, ts, q, o, plot_type, filename=None):

        # TODO: Fix array computation

        poa = self.apm.compute_plot_options_array([ts[0]], [q], [o], plot_type)
        animation = GroupAnimation(ts, poa, self.apm)
        animation.animate()
        animation.save(filename)

    def export_frame(self, t, quantity, orientation, plot_type, filename=None):

        po = self.apm.compute_plot_options(t, quantity, orientation, plot_type)
        if isinstance(po, ScatterPlotOptions):
            plot = Scatter2D(po)
        elif isinstance(po, PColorPlotOptions):
            plot = PColorPlot(po)
        elif isinstance(po, RadialPlotOptions):
            plot = RadialPlot(po)
        else:
            raise ValueError("Unknown plot options type")

        plot.save(filename)

    # def export_group_frames(self, ts, qs, os, plot_type, filename=None):
    #
    #     poa = self.apm.compute_plot_options_array(ts, qs, os, plot_type)
    #     frames = GroupFrame(poa)
    #     frames.save(filename)


if __name__ == "__main__":

    from names import n

    ar = ArepoRun.from_directory("./data/singular_WD")
    aa = ArepoAnalyser(analysis_options={"inner_boxsize": 5e9, "slice_wdith": 0.1})
    apm = PlotManager(ar, aa)
    ah = ArepoHelper(apm)

    for t in range(0, 24):
        po = apm.compute_plot_options(t, n.TEMPERATURE, ["x", "y"], "PColor",
                                      explicit_options={"log_cmap": True})
        plot = PColorPlot(po)
        plot.save(f"{t}.png")

    # ah.export_frame(0, n.TEMPERATURE, ["x", "y"], "PColor")
    # ah.export_animation([0, 45], n.DENSITY, ["x", "y"], "PColor")
    # ah.export_group_frames([7], [n.DENSITY, n.MASSES], [["x", "y"], ["x", "z"]], "PColor")
    # ah.export_group_animation([0, 9], [n.MASSES, n.MAGNETICFIELD], [["x", "y"], ["x", "z"]], "PColor")
