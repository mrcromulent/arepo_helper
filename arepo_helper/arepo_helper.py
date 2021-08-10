from plot_options import PColorPlotOptions, RadialPlotOptions, ScatterPlotOptions
from group_animation import GroupAnimation
from scatter_2d import Scatter2D
from radial_plot import RadialPlot
from pcolor_plot import PColorPlot
from run import ArepoRun
from names import n
from utilities import sci
import numpy as np


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


def main():
    ar = ArepoRun.from_directory("/home/pierre/Desktop/output")
    # ar.make_radial_time_series(n.TEMPERATURE)
    # ar.snapshots[1].quick_radial(n.TEMPERATURE)
    boxsize = ar.snapshots[0].get_from_h5(n.BOXSIZE)
    inner_boxsize = 1e10
    # a = [boxsize / 2 - inner_boxsize / 2, boxsize / 2, boxsize / 2]
    # b = [boxsize / 2 + inner_boxsize / 2, boxsize / 2, boxsize / 2]
    # ar.make_radial_time_series(n.DENSITY, a=a, b=b, cyl_rad=0.01*inner_boxsize, export_filename="density_binary.png")
    # ar.snapshots[-1].quick_pcolor_xnuc(n.NUCLEARCOMPOSITION, inner_boxsize=4e9, export_filename="merger.png",
    #                               select_column=1)

    distances = []
    times = []

    for i, s in enumerate(ar.snapshots):
        passives = s.passive()
        t = s.get_from_h5(n.TIME)
        tol = 1.0e-6
        passives[abs(passives) < tol] = 0.0
        pass0 = passives[:, 0]
        pass1 = passives[:, 1]
        mass = s.mass()

        mass_weighted_pos_wd0 = s.coords() * pass0[:, None] * mass[:, None]
        mass_of_wd0 = (pass0 * mass).sum()

        mass_weighted_pos_wd1 = s.coords() * pass1[:, None] * mass[:, None]
        mass_of_wd1 = (pass1 * mass).sum()

        c0 = mass_weighted_pos_wd0.sum(axis=0) / mass_of_wd0
        c1 = mass_weighted_pos_wd1.sum(axis=0) / mass_of_wd1

        dist = np.sqrt((c1[0] - c0[0]) ** 2 + (c1[1] - c0[1]) ** 2 + (c1[2] - c0[2]) ** 2)
        # print(sci(dist))
        distances.append(dist)
        times.append(t)

    import pandas as pd
    window_size = 3

    numbers_series = pd.Series(distances)
    windows = numbers_series.rolling(window_size)
    moving_averages = windows.mean()

    moving_averages_list = moving_averages.tolist()
    without_nans = moving_averages_list[window_size - 1:]

    print(distances)
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    plt.plot(times, distances)
    # plt.plot(times, moving_averages_list)
    plt.show()
    # for i, s in enumerate(ar.snapshots):
    #     t = s.get_from_h5(n.TIME)
    #     if t > 200:
    #         s.quick_pcolor(n.DENSITY, inner_boxsize=inner_boxsize, export_filename=f"merger_{t}.png")
        # s.quick_nuclear_compositions(export_filename=f"{i}.png")
    # ar.energy_balance()


if __name__ == "__main__":
    main()
