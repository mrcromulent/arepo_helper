from utilities import suppress_stdout_stderr, part_fields
from abstract_plot import AbstractPlot
from names import n
from arepo_vis import make_radial
import numpy as np


class RadialPlot(AbstractPlot):
    line = None

    def __init__(self, plot_options, figure=None, ax=None):
        super(RadialPlot, self).__init__(plot_options, figure, ax)

        tc, tq = self.compute_plot_content()
        self.populate_plot(tc, tq)

    def compute_plot_content(self):

        snap = self.po.ar.snapshots[self.po.t]
        coords = snap.get_from_h5(n.COORDINATES)
        quant = snap.get_from_h5(self.po.quantity)

        tc, tq = self.apply_location_cutoffs(coords, quant)
        tc, tq = self.apply_quantity_cutoffs(tc, tq)

        # Take the magnitude of vector quantities
        if tq.ndim > 1:
            tq = np.apply_along_axis(np.linalg.norm, 1, tq)

        return tc, tq

    def radius(self, coords, center):

        return np.sqrt(
            (coords[:, 0] - center[0]) ** 2 +
            (coords[:, 1] - center[1]) ** 2 +
            (coords[:, 2] - center[2]) ** 2)

    def calc_radial_profile(self, coords, quant):

        center_x = np.average(self.po.xlim)
        center_y = np.average(self.po.ylim)
        center_z = np.average(self.po.zlim)

        min_bs = np.min([self.po.xlim[1] - self.po.xlim[0],
                         self.po.ylim[1] - self.po.ylim[0],
                         self.po.zlim[1] - self.po.zlim[0]])

        a = [center_x, center_y, center_z]
        b = [self.po.xlim[1], center_y, center_z]
        radius = 0.01 * min_bs

        with suppress_stdout_stderr():
            p = make_radial(coords.astype('float64'), quant.astype('float64'),
                            a, b,
                            radius,
                            self.po.nshells)

        return p

    def populate_plot(self, coords, quant):

        p = self.calc_radial_profile(coords, quant)

        if self.po.logscale:
            self.line = self.ax.semilogy(p[1, :], p[0, :], self.po.color)
        else:
            self.line = self.ax.plot(p[1, :], p[0, :], self.po.color)

        units = part_fields[self.po.quantity]["Units"]

        self.set_labels(xlabel=None, ylabel=f"{self.po.quantity} slice [${units}$]")
        self.set_title()
        self.ax.set_ylim(self.po.cbar_lims)

    def update_plot(self, coords, quant):

        super(RadialPlot, self).update_plot()

        p = self.calc_radial_profile(coords, quant)
        line = self.line[0]

        # line.set_data(p[1, :], p[0, :])

        line.set_xdata(p[1, :])
        line.set_ydata(p[0, :])

        self.set_title()
        # self.set_lims(xlim=np.amax(p[1, :]), ylim=np.amax(p[0, :]))
        # self.set_scmp_lims()

    def update_plot_from_plot_options(self, plot_options):

        super(RadialPlot, self).update_plot_from_plot_options()

        self.po = plot_options
        tc, tq = self.compute_plot_content()
        self.update_plot(tc, tq)
