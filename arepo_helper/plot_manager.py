from plot_options import PlotOptions
from names import n, ArepoHeader
from plot_options import mapping
import utilities
import numpy as np


class PlotManager(object):

    def __init__(self, arepo_run, analyser):

        self.ar = arepo_run
        self.aa = analyser

    def compute_plot_options(self, t, quantity, orientation, plot_type, explicit_options=None):

        xlim, ylim, zlim = self.compute_limits()
        xlabel, ylabel = self.compute_labels(orientation)
        title = self.get_title(t, quantity)
        cbar_lims = self.compute_cbar_lims(quantity)

        if explicit_options is not None and "log_cmap" in explicit_options:
            log_cmap = explicit_options["log_cmap"]
        else:
            log_cmap = False

        cmap, norm, scmp = utilities.get_cmap(quantity, cbar_lims, log_cmap=log_cmap)

        plot_options = {"t": t,
                        "ar": self.ar,
                        "aa": self.aa,
                        "xlim": xlim,
                        "ylim": ylim,
                        "zlim": zlim,
                        "title": title,
                        "xlabel": xlabel,
                        "ylabel": ylabel,
                        "quantity": quantity,
                        "orientation": orientation,
                        "cbar_lims": cbar_lims,
                        "scmp": scmp,
                        "cmap": cmap,
                        "norm": norm,
                        "explicit_options": explicit_options}

        if explicit_options is not None:
            for opt in explicit_options.keys():
                plot_options[opt] = explicit_options[opt]

        return mapping[plot_type](plot_options)

    def compute_plot_options_array(self, ts, qs, os, plot_type, explicit_options=None):
        """
        :param ts: list of t values
        :param qs: list of quantities
        :param os: list of orientations
        :param plot_type: string of plot type
        :param explicit_options: dictionary of values to overwrite
        :return: np.array of ArepoPlotOptions objects
        """

        if len(ts) > 1:
            raise NotImplementedError

        t   = ts[0]
        poa = np.empty((len(qs), len(os)), dtype=PlotOptions)

        for i, q in enumerate(qs):
            for j, o in enumerate(os):
                poa[i, j] = self.compute_plot_options(t, q, o, plot_type, explicit_options=explicit_options)

        return poa

    def get_title(self, t, quantity):

        units = utilities.part_fields[quantity]["Units"]
        # time = round(self.ar.snapshots[t].get_from_h5(ArepoHeader.TIME), 2)
        return f"{quantity} time evolution, t = {t} sec. [${units}$]"
        # return f"{quantity} time evolution, t = {time} sec. [${units}$]"

    def compute_labels(self, orientation):
        return f"{orientation[0]} [$cm$]", f"{orientation[1]} [$cm$]"

    def compute_limits(self):

        boxsize = self.ar.run_header[n.BOXSIZE]
        inner_bs = self.find_inner_boxsize()
        sw = self.aa.slice_width

        xlim = [0.5 * boxsize - 0.5 * inner_bs, 0.5 * boxsize + 0.5 * inner_bs]
        ylim = [0.5 * boxsize - 0.5 * inner_bs, 0.5 * boxsize + 0.5 * inner_bs]
        zlim = [0.5 * boxsize - 0.5 * sw * inner_bs, 0.5 * boxsize + 0.5 * sw * inner_bs]

        return xlim, ylim, zlim

    def find_inner_boxsize(self):

        if self.aa.inner_boxsize is None:
            return self.ar.run_header[n.BOXSIZE]
        else:
            return self.aa.inner_boxsize

    def compute_cbar_lims(self, quantity):

        if quantity in utilities.plot_quantities:
            return [np.min(self.ar.minima[quantity]), np.max(self.ar.maxima[quantity])]
        else:
            raise NotImplemented("Not implemented for other quantities.")

    def __str__(self):

        return f"{self.__class__.__name__} \n with ArepoRun: {str(self.ar)} \n and ArepoAnalyser: {str(self.aa)}"
