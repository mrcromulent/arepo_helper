from utilities import SuppressStdout, part_fields, Coordinates as Coords
from mpl_toolkits.axes_grid1 import make_axes_locatable
from arepo_vis import make_radial, make_pcolor
import matplotlib.pyplot as plt
from names import n
import numpy as np
import utilities
import analysis
import run
import os


class PlotOptions(object):
    """Abstract class for plotting options."""

    plot_type = None
    t_idx = 0
    ar = None
    aa = None
    quantity = None
    orientation = None
    xlim = [None, None]
    ylim = [None, None]
    zlim = [None, None]
    title = None
    xlabel = None
    ylabel = None

    def __init__(self, ar, aa, *args, **kwargs):
        """Constructor for abstract plotting options class

        :param ar: Arepo run
        :type ar: run.ArepoRun
        :param aa: Arepo analyser
        :type aa: analysis.ArepoAnalyser
        :param args: Explicit plotting options (e.g. title)
        :param kwargs: Explicit plotting options (e.g. title)

        :rtype: PlotOptions
        :returns: PlotOptions object
        """
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
        """Returns the inner boxsize of the self.ar space"""

        if self.aa.inner_boxsize is None:
            return self.ar.run_header[n.BOXSIZE]
        else:
            return self.aa.inner_boxsize

    def compute_title(self):
        """Computes the title"""
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
        """Abstract constructor of plot options for Scatter2D class

        :param ar: Arepo run
        :type ar: run.ArepoRun
        :param aa: Arepo Analyser
        :type aa: analysis.ArepoAnalyser
        :param args: Explicit options (e.g. title)
        :param kwargs: Explicit options (e.g. title)

        :rtype: ScatterPlotOptions
        :returns: ScatterPlotOptions object
        """

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
        """Computes colorbar limits"""
        if self.quantity in utilities.plot_quantities:
            self.cbar_lims = [np.min(self.ar.minima[self.quantity]), np.max(self.ar.maxima[self.quantity])]
        else:
            raise NotImplemented("Not implemented for other quantities.")

    def compute_limits(self):
        """Computes limits in x,y,z space"""

        boxsize = self.ar.run_header[n.BOXSIZE]
        inner_bs = self.find_inner_boxsize()
        sw = self.slice_width

        self.xlim = [0.5 * boxsize - 0.5 * inner_bs, 0.5 * boxsize + 0.5 * inner_bs]
        self.ylim = [0.5 * boxsize - 0.5 * inner_bs, 0.5 * boxsize + 0.5 * inner_bs]
        self.zlim = [0.5 * boxsize - 0.5 * sw * inner_bs, 0.5 * boxsize + 0.5 * sw * inner_bs]

    def compute_labels(self):
        """Compute x/ylabels"""
        self.xlabel = f"{self.orientation[0]} [$cm$]"
        self.ylabel = f"{self.orientation[1]} [$cm$]"

    def compute_colormap(self):
        """Creates colorbar based on the cbar_lims"""
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
        """Plotting options for RadialPlot

        :param ar: Arepo run
        :type ar: run.ArepoRun
        :param aa: Arepo Analyser
        :type aa: analysis.ArepoAnalyser
        :param args: Explicit plot options (e.g. title)
        :param kwargs: Explicit plot options (e.g. title)

        :rtype: RadialPlotOptions
        :return: RadialPlotOptions object
        """

        super(RadialPlotOptions, self).__init__(ar, aa, *args, **kwargs)

        self.compute_title()
        self.compute_labels()
        self.compute_ab_radius()

        for key in aa.analysis_options:
            if aa.analysis_options[key] is not None:
                setattr(self, key, aa.analysis_options[key])

    def compute_ab_radius(self):
        """Computes the a & b points of the radial cylinder, as well as its radius"""

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
        """Computes the x/ylabels of the plot"""
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
        """Plotting options for PColorPlot

        :param ar: Arepo Run
        :type ar: run.ArepoRun
        :param aa: Arepo Analyser
        :type aa: analysis.ArepoAnalyser
        :param args: Explicit plot options (e.g. title)
        :param kwargs: Explicit plot options (e.g. title)

        :rtype: PColorPlotOptions
        :returns: PColorPlotOptions object
        """

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
        """Computes the colormap with limits cbar_lims"""
        self.cmap, self.norm, self.scmp = utilities.get_cmap(self.quantity, self.cbar_lims, log_cmap=self.log_cmap)

    def compute_cbar_lims(self):
        """Computes the cbar_lims of the colorbar"""
        if self.quantity in utilities.plot_quantities:
            self.cbar_lims = [np.min(self.ar.minima[self.quantity]), np.max(self.ar.maxima[self.quantity])]
        else:
            raise NotImplementedError(f"Not implemented for {self.quantity}")

    def compute_limits(self):
        """Computes the x,y,z limits of the self.ar space"""

        boxsize = self.ar.snapshots[0].get_from_h5(n.BOXSIZE)
        inner_bs = self.find_inner_boxsize()

        self.xlim = [0.5 * boxsize - 0.5 * inner_bs, 0.5 * boxsize + 0.5 * inner_bs]
        self.ylim = [0.5 * boxsize - 0.5 * inner_bs, 0.5 * boxsize + 0.5 * inner_bs]
        self.zlim = [0.5 * boxsize - 0.5 * inner_bs, 0.5 * boxsize + 0.5 * inner_bs]

    def compute_labels(self):
        """Computes the x/ylabels"""
        self.xlabel = f"{self.orientation[0]} [$cm$]"
        self.ylabel = f"{self.orientation[1]} [$cm$]"


class AbstractPlot(object):
    """Abstract class to for Arepo type plots
    """

    cb = None

    def __init__(self, plot_options, figure=None, ax=None):
        """Constructor for abstract plot class

        :param plot_options: A plot options object which contains information about what appears on the plot
        :type plot_options: Any
        :param figure: Optional matplotlib figure to plot results on
        :type figure: matplotlib.figure.Figure
        :param ax: Optional matplotlib axes to plot results on
        :type ax: matplotlib.axes._subplots.AxesSubplot

        :rtype: AbstractPlot
        :returns: AbstractPlot object
        """

        self.po = plot_options

        if figure is None or ax is None:
            self.fig, self.ax = plt.subplots()
        else:
            self.fig, self.ax = figure, ax

    def set_title(self, new_title):
        """Sets the title"""
        self.ax.set_title(new_title)

    def set_labels(self, xlabel, ylabel):
        """Sets the x/ylabels

        :param xlabel:
        :type xlabel: str
        :param ylabel:
        :type ylabel: str
        """
        self.ax.set(xlabel=xlabel, ylabel=ylabel)

    def set_lims(self, xlim, ylim):
        """Sets the limits

        :param xlim:
        :type xlim: list
        :param ylim:
        :type ylim: list
        """
        self.ax.set(xlim=xlim, ylim=ylim)

    def apply_location_cutoffs(self, coords, quant):
        """Applies location cutoffs in the space defined by an arepo run

        :param coords: Coordinate array
        :type coords: np.ndarray
        :param quant: Quantity array
        :type quant: np.ndarray
        :rtype: (np.ndarray, np.ndarray)
        :return: Trimmed coordinates and trimmed quantity
        """
        xc, yc, zc = Coords.coordinates(self.po.orientation)

        mask = np.multiply(coords[:, zc] > self.po.zlim[0], coords[:, zc] < self.po.zlim[1])
        trimmed_coords = coords[mask, :]
        trimmed_quant = quant[mask]

        xmask = np.multiply(trimmed_coords[:, xc] > self.po.xlim[0], trimmed_coords[:, xc] < self.po.xlim[1])
        ymask = np.multiply(trimmed_coords[:, yc] > self.po.ylim[0], trimmed_coords[:, yc] < self.po.ylim[1])
        mask = np.multiply(xmask, ymask)

        return trimmed_coords[mask, :], trimmed_quant[mask]

    def apply_quantity_cutoffs(self, coords, quant):
        """Removes all but the points which fall in the ranges specified by the cutoff_table

        :param coords: Coordinate array
        :type coords: np.ndarray
        :param quant: Quantity array
        :type quant: np.ndarray
        :rtype: (np.ndarray, np.ndarray)
        :return: Trimmed coordinates and trimmed quantity
        """

        if self.po.quantity in self.po.aa.cutoff_table.keys():
            cutoffs = self.po.aa.cutoff_table[self.po.quantity]
            mask = np.multiply(quant > cutoffs[0], quant < cutoffs[1])
            coords = coords[mask, :]
            quant = quant[mask]

        return coords, quant

    def save(self, filename=None):
        """Saves the plot to the location specified by filename

        :param filename: None or absolute filename or file stem
        :type filename: str
        """

        if filename is None:
            fn_stem = f"{type(self).__name__}_{self.po.t_idx}_{self.po.quantity}_{self.po.orientation}.png"
            dir_name = os.path.dirname(self.po.ar.snapshots[0].filename)
            filename = os.path.join(dir_name, fn_stem)
        elif not os.path.isabs(filename):
            dir_name = os.path.dirname(self.po.ar.snapshots[0].filename)
            filename = os.path.join(dir_name, filename)

        self.fig.savefig(filename, dpi=300)

    def populate_plot(self):
        """Populates the plot

        :return:
        """
        pass

    def update_plot(self):
        pass

    def compute_plot_content(self):
        """Computes the content of the plot using location and quantity cutoffs

        :return: Trimmed coordinates and quantities
        :rtype: (np.ndarray, np.ndarray)
        """
        pass

    def update_plot_from_plot_options(self):
        pass


class RadialPlot(AbstractPlot):

    def __init__(self, plot_options, figure=None, ax=None, do_not_plot=False):
        """Radial Plot

        :param plot_options: Plot options
        :type plot_options: RadialPlotOptions
        :param figure: Optional matplotlib figure to plot results on
        :type figure: matplotlib.figure.Figure
        :param ax: Optional matplotlib axes to plot results on
        :type ax: matplotlib.axes._subplots.AxesSubplot
        :param do_not_plot: Stops the plot from being plotted
        :type do_not_plot: bool

        :rtype: RadialPlot
        :return: RadialPlot object
        """

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
        """Calculates the content of the cylinder

        :param coords: Coordinates
        :type coords: np.ndarray
        :param quant: Quantity
        :type quant: np.ndarray
        :return: p
        :rtype: np.ndarray
        """

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
        """Updates the plot from plot options

        :param plot_options: New plot options
        :type plot_options: RadialPlotOptions
        """

        super(RadialPlot, self).update_plot_from_plot_options()

        self.po = plot_options
        tc, tq = self.compute_plot_content()
        self.update_plot(tc, tq)


class Scatter2D(AbstractPlot):
    """Class for scatter plots
    """

    scatter = None
    cb = None

    def __init__(self, plot_options, figure=None, ax=None):
        """Constructor

        :param plot_options: Plot options for scatter
        :type plot_options: ScatterPlotOptions
        :param figure: Optional matplotlib figure to plot results on
        :type figure: matplotlib.figure.Figure
        :param ax: Optional matplotlib axes to plot results on
        :type ax: matplotlib.axes._subplots.AxesSubplot

        :rtype: Scatter2D
        :returns: Scatter2D object
        """

        super(Scatter2D, self).__init__(plot_options, figure, ax)

        tc, tq = self.compute_plot_content()
        self.populate_plot(tc, tq)

    def set_a(self, quant):
        """Sets the plot frequency. i.e. every ath point is plotted"""

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

    def set_scmp_lims(self, lims):
        """Sets limits on the scalar mappable"""
        self.po.scmp.set_clim(lims)

    def set_colorbar(self):
        """Creates the colorbar"""
        if self.po.include_colorbar:
            divider = make_axes_locatable(self.ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            self.cb = self.fig.colorbar(self.po.scmp, cax=cax)
            self.set_scmp_lims(self.po.cbar_lims)
            self.cb.set_alpha(1)
            self.cb.draw_all()

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
        """Updates the plot from plot options

        :param plot_options: New plot options
        :type plot_options: ScatterPlotOptions
        """

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
    """Class to create PColors
    """

    pcolor = None

    def __init__(self, plot_options, figure=None, ax=None, do_not_plot=False):
        """PColor Plot constructor

        :param plot_options: Plot options
        :type plot_options: PColorPlotOptions
        :param figure: Optional matplotlib figure to plot results on
        :type figure: matplotlib.figure.Figure
        :param ax: Optional matplotlib axes to plot results on
        :type ax: matplotlib.axes._subplots.AxesSubplot
        :param do_not_plot: Stops the plot from being plotted
        :type do_not_plot: bool

        :rtype: PColorPlot
        :return: PColorPlot object
        """

        super(PColorPlot, self).__init__(plot_options, figure, ax)

        self.x = None
        self.y = None
        self.data = None

        self.compute_plot_content()
        if not do_not_plot:
            self.populate_plot()

    def calc_a_slice(self, coords, quant):
        """Computes the content of the pcolor

        :param coords: Coordinates
        :type coords: np.ndarray
        :param quant: Quantity
        :type quant: np.ndarray
        :return: x, y and data dict containing values under the key "grid"
        :rtype: dict
        """

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

    def set_scmp_lims(self, lims):
        """Sets limits on the scalar mappable"""
        self.po.scmp.set_clim(lims)

    def set_colorbar(self):
        """Creates the colorbar"""
        if self.po.include_colorbar:
            divider = make_axes_locatable(self.ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            self.cb = self.fig.colorbar(self.po.scmp, cax=cax)
            self.set_scmp_lims(self.po.cbar_lims)
            self.cb.set_alpha(1)
            self.cb.draw_all()

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
        """Updates the plot from plot options

        :param plot_options: New plot options
        :type plot_options: PColorPlotOptions
        """

        super(PColorPlot, self).update_plot_from_plot_options()

        self.po = plot_options
        tc, tq = self.compute_plot_content()
        self.update_plot(tc, tq)


mapping = {ScatterPlotOptions.plot_type: [ScatterPlotOptions, Scatter2D],
           PColorPlotOptions.plot_type: [PColorPlotOptions, PColorPlot],
           RadialPlotOptions.plot_type: [RadialPlotOptions, RadialPlot]}


def get_plotter_func(plot_name):
    """Gets the corresponding Plot class from the string name

    :param plot_name: string name of plot class
    :type plot_name: str
    :return: Plot class
    :rtype: Any
    """
    return mapping[plot_name][1]


def get_plotter_func_from_plot_options(po):
    """Gets the corresponding Plot class from an instance of plot options.

    :param po: plot_options
    :type po: Any
    :return: Plot class
    :rtype: Any
    """
    for plot_type in mapping:
        if isinstance(po, mapping[plot_type][0]):
            return mapping[plot_type][1]

    raise ValueError("Unknown plot options type")


def get_plot_options_func(plot_name):
    """Gets the corresponding Plot Options class from the plot name.

    :param plot_name: string name of plot class
    :type plot_name: str
    :return: Plot class
    :rtype: Any
    """
    return mapping[plot_name][0]
