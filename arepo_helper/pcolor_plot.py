from utilities import suppress_stdout_stderr
from abstract_plot import AbstractPlot
from utilities import Coordinates as C
from arepo_vis import make_pcolor
from names import n
import numpy as np


class PColorPlot(AbstractPlot):
    pcolor = None

    def __init__(self, plot_options, figure=None, ax=None):
        super(PColorPlot, self).__init__(plot_options, figure, ax)

        tc, tq = self.compute_plot_content()
        self.populate_plot(tc, tq)

    def calc_a_slice(self, coords, quant):
        res = self.po.resolution
        boxsize_x = self.po.xlim[1] - self.po.xlim[0]
        boxsize_y = self.po.ylim[1] - self.po.ylim[0]
        numthreads = self.po.numthreads

        resolutions = [res, res]
        boxsizes = [boxsize_x, boxsize_y]
        centers = [np.average(self.po.xlim),
                   np.average(self.po.ylim),
                   np.average(self.po.zlim)]
        xc, yc, _ = C.coordinates(self.po.orientation)
        axes = [xc, yc]

        with suppress_stdout_stderr():
            data = make_pcolor(coords.astype('float64'), quant.astype('float64'),
                               axes,
                               boxsizes,
                               resolutions,
                               centers,
                               include_neighbours_in_output=1,
                               numthreads=numthreads)

        x = np.arange(res + 1, dtype="float64") / res * boxsize_x - 0.5 * boxsize_x + centers[0]
        y = np.arange(res + 1, dtype="float64") / res * boxsize_y - 0.5 * boxsize_y + centers[1]

        return x, y, data

    def populate_plot(self, coords, points):
        super(PColorPlot, self).populate_plot()

        x, y, data = self.calc_a_slice(coords, points)
        self.pcolor = self.ax.pcolormesh(x, y, np.transpose(data["grid"]), shading='flat')
        self.set_labels()
        self.set_title()
        self.set_lims()
        self.set_colorbar()

        return self.pcolor

    def compute_plot_content(self):
        super(PColorPlot, self).compute_plot_content()

        snap = self.po.ar.snapshots[self.po.t]
        coords = snap.get_from_h5(n.COORDINATES)
        quant = snap.get_from_h5(self.po.quantity)

        tc, tq = self.apply_location_cutoffs(coords, quant)
        tc, tq = self.apply_quantity_cutoffs(tc, tq)

        # Take the magnitude of vector quantities
        if tq.ndim > 1:
            tq = np.apply_along_axis(np.linalg.norm, 1, tq)
            print(tq)

        return tc, tq

    def update_plot(self, coords, quant):
        super(PColorPlot, self).update_plot()

        x, y, data = self.calc_a_slice(coords, quant)
        z = np.transpose(data["grid"])

        self.pcolor.set_array(z.ravel())
        self.set_title()
        self.set_scmp_lims()

    def update_plot_from_plot_options(self, plot_options):
        super(PColorPlot, self).update_plot_from_plot_options()

        self.po = plot_options
        tc, tq = self.compute_plot_content()
        self.update_plot(tc, tq)
