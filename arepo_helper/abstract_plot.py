from utilities import SuppressStdout, part_fields, Coordinates as Coords
from mpl_toolkits.axes_grid1 import make_axes_locatable
from arepo_vis import make_radial, make_pcolor
import matplotlib.pyplot as plt
from names import n
import numpy as np
import utilities
import os


class PlotOptions(object):

    plot_type = None
    t_idx = 0
    ar = None
    aa = None
    quantity = None
    orientation = None
    xlim = [None, None]
    ylim = [None, None]
    title = None
    xlabel = None
    ylabel = None

    def __init__(self, ar, aa, *args, **kwargs):

        self.ar = ar
        self.aa = aa

        for key in aa.analysis_options:
            setattr(self, key, aa.analysis_options[key])

        for dictionary in args:
            for key in dictionary:
                setattr(self, key, dictionary[key])
        for key in kwargs:
            setattr(self, key, kwargs[key])

    def find_inner_boxsize(self):

        if self.aa.inner_boxsize is None:
            return self.ar.run_header[n.BOXSIZE]
        else:
            return self.aa.inner_boxsize

    def compute_title(self):
        units = part_fields[self.quantity]["Units"]
        t = round(self.ar.snapshots[self.t_idx].get_from_h5(n.TIME), 2)
        self.title = f"{self.quantity} time evolution, t = {t} sec. [${units}$]"


class ScatterPlotOptions(PlotOptions):

    plot_type           = "Scatter"
    psize               = 0.01
    a                   = 2
    slice_width         = 0.005
    cmap                = None
    scmp                = None
    norm                = None
    include_colorbar    = True
    cbar_lims           = [None, None]
    zlim                = [None, None]
    log_cmap            = True

    def __init__(self, ar, aa, *args, **kwargs):
        super(ScatterPlotOptions, self).__init__(ar, aa, *args, **kwargs)

        self.compute_title()
        self.compute_limits()
        self.compute_labels()
        self.compute_cbar_lims()
        self.compute_colormap()

        for key in aa.analysis_options:
            if aa.analysis_options[key] is not None:
                setattr(self, key, aa.analysis_options[key])

    def compute_cbar_lims(self):
        if self.quantity in utilities.plot_quantities:
            self.cbar_lims = [np.min(self.ar.minima[self.quantity]), np.max(self.ar.maxima[self.quantity])]
        else:
            raise NotImplemented("Not implemented for other quantities.")

    def compute_limits(self):

        boxsize = self.ar.run_header[n.BOXSIZE]
        inner_bs = self.find_inner_boxsize()
        sw = self.slice_width

        self.xlim = [0.5 * boxsize - 0.5 * inner_bs, 0.5 * boxsize + 0.5 * inner_bs]
        self.ylim = [0.5 * boxsize - 0.5 * inner_bs, 0.5 * boxsize + 0.5 * inner_bs]
        self.zlim = [0.5 * boxsize - 0.5 * sw * inner_bs, 0.5 * boxsize + 0.5 * sw * inner_bs]

    def compute_labels(self):
        self.xlabel = f"{self.orientation[0]} [$cm$]"
        self.ylabel = f"{self.orientation[1]} [$cm$]"

    def compute_colormap(self):
        self.cmap, self.norm, self.scmp = utilities.get_cmap(self.quantity, self.cbar_lims, log_cmap=self.log_cmap)


class RadialPlotOptions(PlotOptions):

    plot_type           = "Radial"
    logscale            = False
    nshells             = 200
    color               = "k"
    include_colorbar    = False
    a                   = None
    b                   = None
    radius              = None
    zlim                = None

    def __init__(self, ar, aa, *args, **kwargs):
        super(RadialPlotOptions, self).__init__(ar, aa, *args, **kwargs)

        self.compute_title()
        self.compute_labels()
        self.compute_ab_radius()

        for key in aa.analysis_options:
            if aa.analysis_options[key] is not None:
                setattr(self, key, aa.analysis_options[key])

    def compute_ab_radius(self):

        inner_bs = self.find_inner_boxsize()
        boxsize = self.ar.run_header[n.BOXSIZE]
        self.radius = 0.01 * inner_bs

        self.xlim = [0.5 * boxsize - 0.5 * inner_bs, 0.5 * boxsize + 0.5 * inner_bs]
        self.ylim = [0.5 * boxsize - 0.5 * inner_bs, 0.5 * boxsize + 0.5 * inner_bs]
        self.zlim = [0.5 * boxsize - 0.5 * inner_bs, 0.5 * boxsize + 0.5 * inner_bs]

        if self.a is None or self.b is None:
            center_x = np.average(self.xlim)
            center_y = np.average(self.ylim)
            center_z = np.average(self.zlim)

            #
            self.a = [center_x, center_y, center_z]
            self.b = [self.xlim[1], center_y, center_z]

    def compute_labels(self):
        self.xlabel = f"r [$cm$]"
        self.ylabel = f"{self.quantity}"


class PColorPlotOptions(PlotOptions):

    plot_type           = "PColor"
    resolution          = 1024
    numthreads          = 1
    include_colorbar    = True
    cmap                = None
    scmp                = None
    norm                = None
    cbar_lims           = [None, None]
    zlim                = [None, None]
    log_cmap            = True
    select_column       = None

    def __init__(self, ar, aa, *args, **kwargs):
        super(PColorPlotOptions, self).__init__(ar, aa, *args, **kwargs)

        self.compute_title()
        self.compute_limits()
        self.compute_labels()
        self.compute_cbar_lims()
        self.compute_colormap()

        for key in aa.analysis_options:
            if aa.analysis_options[key] is not None:
                setattr(self, key, aa.analysis_options[key])

    def compute_colormap(self):
        self.cmap, self.norm, self.scmp = utilities.get_cmap(self.quantity, self.cbar_lims, log_cmap=self.log_cmap)

    def compute_cbar_lims(self):
        if self.quantity in utilities.plot_quantities:
            self.cbar_lims = [np.min(self.ar.minima[self.quantity]), np.max(self.ar.maxima[self.quantity])]
        else:
            raise NotImplementedError(f"Not implemented for {self.quantity}")

    def compute_limits(self):

        boxsize = self.ar.run_header[n.BOXSIZE]
        inner_bs = self.find_inner_boxsize()

        self.xlim = [0.5 * boxsize - 0.5 * inner_bs, 0.5 * boxsize + 0.5 * inner_bs]
        self.ylim = [0.5 * boxsize - 0.5 * inner_bs, 0.5 * boxsize + 0.5 * inner_bs]
        self.zlim = [0.5 * boxsize - 0.5 * inner_bs, 0.5 * boxsize + 0.5 * inner_bs]

    def compute_labels(self):
        self.xlabel = f"{self.orientation[0]} [$cm$]"
        self.ylabel = f"{self.orientation[1]} [$cm$]"


class AbstractPlot(object):

    cb = None

    def __init__(self, plot_options, figure=None, ax=None):
        self.po = plot_options

        if figure is None:
            self.fig, self.ax = plt.subplots()
        else:
            self.fig, self.ax = figure, ax

    def set_title(self, new_title):
        self.ax.set_title(new_title)

    def set_labels(self, xlabel, ylabel):
        self.ax.set(xlabel=xlabel, ylabel=ylabel)

    def set_lims(self, xlim, ylim):
        self.ax.set(xlim=xlim, ylim=ylim)

    def set_scmp_lims(self, lims):
        self.po.scmp.set_clim(lims)

    def set_colorbar(self):
        if self.po.include_colorbar:
            divider = make_axes_locatable(self.ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            self.cb = self.fig.colorbar(self.po.scmp, cax=cax)
            self.set_scmp_lims(self.po.cbar_lims)
            self.cb.set_alpha(1)
            self.cb.draw_all()

    def apply_location_cutoffs(self, coords, quant):
        xc, yc, zc = Coords.coordinates(self.po.orientation)

        mask = np.multiply(coords[:, zc] > self.po.zlim[0], coords[:, zc] < self.po.zlim[1])
        trimmed_coords = coords[mask, :]
        trimmed_quant = quant[mask]

        xmask = np.multiply(trimmed_coords[:, xc] > self.po.xlim[0], trimmed_coords[:, xc] < self.po.xlim[1])
        ymask = np.multiply(trimmed_coords[:, yc] > self.po.ylim[0], trimmed_coords[:, yc] < self.po.ylim[1])
        mask = np.multiply(xmask, ymask)

        return trimmed_coords[mask, :], trimmed_quant[mask]

    def apply_quantity_cutoffs(self, trimmed_coords, trimmed_quant):

        if self.po.quantity in self.po.aa.cutoff_table.keys():
            cutoffs = self.po.aa.cutoff_table[self.po.quantity]
            mask = np.multiply(trimmed_quant > cutoffs[0], trimmed_quant < cutoffs[1])

            trimmed_coords = trimmed_coords[mask, :]
            trimmed_quant = trimmed_quant[mask]

        return trimmed_coords, trimmed_quant

    def save(self, filename=None):

        if filename is None:
            fn_stem = f"{type(self).__name__}_{self.po.t_idx}_{self.po.quantity}_{self.po.orientation}.png"
            dir_name = os.path.dirname(self.po.ar.snapshots[0].filename)
            filename = os.path.join(dir_name, fn_stem)
        elif not os.path.isabs(filename):
            dir_name = os.path.dirname(self.po.ar.snapshots[0].filename)
            filename = os.path.join(dir_name, filename)

        self.fig.savefig(filename, dpi=300)

    def populate_plot(self):
        pass

    def update_plot(self):
        pass

    def compute_plot_content(self):
        pass

    def update_plot_from_plot_options(self):
        pass


class RadialPlot(AbstractPlot):

    def __init__(self, plot_options, figure=None, ax=None, do_not_plot=False):
        super(RadialPlot, self).__init__(plot_options, figure, ax)

        self.line    = None
        self.x       = None
        self.y       = None

        self.compute_plot_content()
        if not do_not_plot:
            self.populate_plot()

    def compute_plot_content(self):

        snap = self.po.ar.snapshots[self.po.t_idx]
        coords = snap.get_from_h5(n.COORDINATES)
        quant = snap.get_from_h5(self.po.quantity)

        tc, tq = self.apply_location_cutoffs(coords, quant)
        tc, tq = self.apply_quantity_cutoffs(tc, tq)

        # Take the magnitude of vector quantities
        if tq.ndim > 1:
            if self.po.aa.select_column is not None:
                tq = tq[:, self.po.aa.select_column]
            else:
                tq = np.apply_along_axis(np.linalg.norm, 1, tq)

        p = self.calc_radial_profile(tc, tq)
        self.x, self.y = p[1, :], p[0, :]

        return tc, tq

    def calc_radial_profile(self, coords, quant):

        with SuppressStdout():
            p = make_radial(coords.astype("float64"), quant.astype("float64"),
                            self.po.a, self.po.b,
                            self.po.radius,
                            self.po.nshells)

        return p

    def populate_plot(self):

        if self.po.logscale:
            self.line = self.ax.semilogy(self.x, self.y, color=self.po.color)
        else:
            self.line = self.ax.plot(self.x, self.y, color=self.po.color)

        self.set_labels(self.po.xlabel, self.po.ylabel)
        self.set_title(self.po.title)

    def update_plot(self, coords, quant):

        super(RadialPlot, self).update_plot()

        p = self.calc_radial_profile(coords, quant)
        self.line[0].set_xdata(p[1, :])
        self.line[0].set_ydata(p[0, :])
        self.set_title(self.po.title)

    def update_plot_from_plot_options(self, plot_options):

        super(RadialPlot, self).update_plot_from_plot_options()

        self.po = plot_options
        tc, tq = self.compute_plot_content()
        self.update_plot(tc, tq)


class Scatter2D(AbstractPlot):
    scatter = None
    cb = None

    def __init__(self, plot_options, figure=None, ax=None):
        super(Scatter2D, self).__init__(plot_options, figure, ax)

        tc, tq = self.compute_plot_content()
        self.populate_plot(tc, tq)

    def set_a(self, quant):
        num_points = np.size(quant, axis=0)
        self.po.a = int(np.ceil(num_points / 100_000))

    def populate_plot(self, coords, quant):

        super(Scatter2D, self).populate_plot()

        xc, yc, _ = Coords.coordinates(self.po.orientation)
        x = coords[:, xc]
        y = coords[:, yc]

        self.scatter = self.ax.scatter(x, y,
                                       c=quant,
                                       s=self.po.psize,
                                       cmap=self.po.cmap,
                                       norm=self.po.norm)

        self.set_labels(self.po.xlabel, self.po.ylabel)
        self.set_title(self.po.title)
        self.set_lims(self.po.xlim, self.po.ylim)
        self.set_colorbar()
        self.set_a(quant)

        return self.scatter

    def update_plot(self, coords, quant):

        super(Scatter2D, self).update_plot()

        xc, yc, zc = Coords.coordinates(self.po.orientation)
        coords = np.delete(coords, zc, axis=1)

        self.set_a(quant)
        self.scatter.set_offsets(coords[::self.po.a])
        self.scatter.set_array(quant[::self.po.a])

        self.set_title(self.po.title)
        self.set_scmp_lims(self.po.cbar_lims)

    def update_plot_from_plot_options(self, plot_options):

        super(Scatter2D, self).update_plot_from_plot_options()

        self.po = plot_options
        tc, tq = self.compute_plot_content()
        self.update_plot(tc, tq)

    def compute_plot_content(self):

        super(Scatter2D, self).compute_plot_content()

        snap = self.po.ar.snapshots[self.po.t_idx]
        coords = snap.get_from_h5(n.COORDINATES)
        quant = snap.get_from_h5(self.po.quantity)

        tc, tq = self.apply_location_cutoffs(coords, quant)
        tc, tq = self.apply_quantity_cutoffs(tc, tq)

        # Take the magnitude of vector quantities
        if tq.ndim > 1:
            tq = np.apply_along_axis(np.linalg.norm, 1, tq)

        return tc, tq


class PColorPlot(AbstractPlot):
    pcolor = None

    def __init__(self, plot_options, figure=None, ax=None, do_not_plot=False):
        super(PColorPlot, self).__init__(plot_options, figure, ax)

        self.x = None
        self.y = None
        self.data = None

        self.compute_plot_content()
        if not do_not_plot:
            self.populate_plot()

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
        xc, yc, _ = Coords.coordinates(self.po.orientation)
        axes = [xc, yc]

        with SuppressStdout():
            data = make_pcolor(coords.astype("float64"), quant.astype("float64"),
                               axes,
                               boxsizes,
                               resolutions,
                               centers,
                               include_neighbours_in_output=1,
                               numthreads=numthreads)

        x = np.arange(res + 1, dtype="float64") / res * boxsize_x - 0.5 * boxsize_x + centers[0]
        y = np.arange(res + 1, dtype="float64") / res * boxsize_y - 0.5 * boxsize_y + centers[1]

        return x, y, data

    def populate_plot(self):
        super(PColorPlot, self).populate_plot()

        self.pcolor = self.ax.pcolormesh(self.x, self.y,
                                         np.transpose(self.data["grid"]),
                                         shading="flat", rasterized=True,
                                         norm=self.po.norm, cmap=self.po.cmap)
        self.set_labels(self.po.xlabel, self.po.ylabel)
        self.set_title(self.po.title)
        self.set_lims(self.po.xlim, self.po.ylim)
        self.set_colorbar()

        return self.pcolor

    def compute_plot_content(self):
        super(PColorPlot, self).compute_plot_content()

        snap = self.po.ar.snapshots[self.po.t_idx]
        coords = snap.get_from_h5(n.COORDINATES)
        quant = snap.get_from_h5(self.po.quantity)

        tc, tq = self.apply_location_cutoffs(coords, quant)
        tc, tq = self.apply_quantity_cutoffs(tc, tq)

        # Take the magnitude of vector quantities
        if tq.ndim > 1:
            if self.po.aa.select_column is not None:
                tq = tq[:, self.po.aa.select_column]
            else:
                tq = np.apply_along_axis(np.linalg.norm, 1, tq)

        self.x, self.y, self.data = self.calc_a_slice(tc, tq)

        return tc, tq

    def update_plot(self, coords, quant):
        super(PColorPlot, self).update_plot()

        self.x, self.y, self.data = self.calc_a_slice(coords, quant)
        z = np.transpose(self.data["grid"])

        self.pcolor.set_array(z.ravel())
        self.set_title(self.po.title)
        self.set_scmp_lims(self.po.cbar_lims)

    def update_plot_from_plot_options(self, plot_options):
        super(PColorPlot, self).update_plot_from_plot_options()

        self.po = plot_options
        tc, tq = self.compute_plot_content()
        self.update_plot(tc, tq)


mapping = {ScatterPlotOptions.plot_type: [ScatterPlotOptions, Scatter2D],
           PColorPlotOptions.plot_type: [PColorPlotOptions, PColorPlot],
           RadialPlotOptions.plot_type: [RadialPlotOptions, RadialPlot]}


def get_plotter_func(plot_name):
    return mapping[plot_name][1]


def get_plotter_func_from_plot_options(po):
    for plot_type in mapping:
        if isinstance(po, mapping[plot_type][0]):
            return mapping[plot_type][1]

    raise ValueError("Unknown plot options type")


def get_plot_options_func(plot_name):
    return mapping[plot_name][0]
