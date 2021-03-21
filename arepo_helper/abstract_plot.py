from mpl_toolkits.axes_grid1 import make_axes_locatable
from utilities import Coordinates as C
import matplotlib.pyplot as plt
import numpy as np


class AbstractPlot(object):

    cb = None

    def __init__(self, plot_options, figure=None, ax=None):
        self.po = plot_options

        if figure is None:
            self.fig, self.ax = plt.subplots()
        else:
            self.fig, self.ax = figure, ax

    def set_title(self, new_title=None):
        if new_title is None:
            new_title = self.po.title

        self.ax.set_title(new_title)

    def set_labels(self, xlabel=None, ylabel=None):
        if xlabel is None:
            xlabel = self.po.xlabel
        if ylabel is None:
            ylabel = self.po.ylabel

        self.ax.set(xlabel=xlabel, ylabel=ylabel)

    def set_lims(self, xlim=None, ylim=None):
        if xlim is None:
            xlim = self.po.xlim
        if ylim is None:
            ylim = self.po.ylim

        self.ax.set(xlim=xlim, ylim=ylim)

    def set_scmp_lims(self, lims=None):
        if lims is None:
            lims = self.po.cbar_lims

        self.po.scmp.set_clim(lims)

    def set_colorbar(self):

        if self.po.include_colorbar:
            divider = make_axes_locatable(self.ax)
            cax = divider.append_axes('right', size='5%', pad=0.05)
            self.cb = self.fig.colorbar(self.po.scmp, cax=cax)
            self.set_scmp_lims()
            self.cb.set_alpha(1)
            self.cb.draw_all()

    def apply_location_cutoffs(self, coords, quant):
        po = self.po
        xc, yc, zc = C.coordinates(self.po.orientation)

        mask = np.multiply(coords[:, zc] > po.zlim[0], coords[:, zc] < po.zlim[1])
        trimmed_coords = coords[mask, :]
        trimmed_quant = quant[mask]

        xmask = np.multiply(trimmed_coords[:, xc] > po.xlim[0], trimmed_coords[:, xc] < po.xlim[1])
        ymask = np.multiply(trimmed_coords[:, yc] > po.ylim[0], trimmed_coords[:, yc] < po.ylim[1])
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
            filename = f"./{type(self).__name__}_{self.po.t}_{self.po.quantity}_{self.po.orientation}.png"

        self.fig.savefig(filename, dpi=300)

    def populate_plot(self):
        pass

    def update_plot(self):
        pass

    def compute_plot_content(self):
        pass

    def update_plot_from_plot_options(self):
        pass
