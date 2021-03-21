from abstract_plot import AbstractPlot
from utilities import Coordinates as C
from names import n
import numpy as np


class Scatter2D(AbstractPlot):
    highlighted_scatter = None
    scatter = None

    def __init__(self, plot_options, figure=None, ax=None):
        super(Scatter2D, self).__init__(plot_options, figure, ax)

        tc, tq, hc, hq = self.compute_plot_content()
        self.populate_plot(tc, tq, hc, hq)

    def set_alpha(self, highlight_coords):
        if highlight_coords is not None:
            self.scatter.set_alpha(self.po.alpha_excluded)

    def set_a(self, quant):
        num_points = np.size(quant, axis=0)
        self.po.a = int(np.ceil(num_points / 100_000))

    def add_highlights(self, highlight_coords, highlight_quant):
        if highlight_coords is not None:
            xc, yc, _ = C.coordinates(self.po.orientation)
            x = highlight_coords[:, xc]
            y = highlight_coords[:, yc]

            self.highlighted_scatter = self.ax.scatter(x, y,
                                                       c=highlight_quant,
                                                       s=self.po.psize,
                                                       alpha=self.po.alpha_included,
                                                       cmap=self.po.cmap,
                                                       norm=self.po.norm)

    def apply_highlight_cutoffs(self, trimmed_coords, trimmed_quant):

        if self.po.quantity in self.po.aa.highlights.keys():
            cutoffs = self.po.aa.highlights[self.po.quantity]
            mask = np.multiply(trimmed_quant > cutoffs[0], trimmed_quant < cutoffs[1])

            hc = trimmed_coords[mask, :]
            hq = trimmed_quant[mask]

            trimmed_coords = trimmed_coords[~mask, :]
            trimmed_quant = trimmed_quant[~mask]

            return trimmed_coords, trimmed_quant, hc, hq

        return trimmed_coords, trimmed_quant, None, None

    def populate_plot(self, coords, quant, highlight_coords=None, highlight_quant=None):

        super(Scatter2D, self).populate_plot()

        xc, yc, _ = C.coordinates(self.po.orientation)
        x = coords[:, xc]
        y = coords[:, yc]

        self.scatter = self.ax.scatter(x, y,
                                       c=quant,
                                       s=self.po.psize,
                                       alpha=self.po.alpha_included,
                                       cmap=self.po.cmap,
                                       norm=self.po.norm)

        self.add_highlights(highlight_coords, highlight_quant)
        self.set_alpha(highlight_coords)
        self.set_labels()
        self.set_title()
        self.set_lims()
        self.set_colorbar()
        self.set_a(quant)

        return self.scatter

    def update_plot(self, coords, quant, highlight_coords=None, highlight_quant=None):

        super(Scatter2D, self).update_plot()

        xc, yc, zc = C.coordinates(self.po.orientation)
        coords = np.delete(coords, zc, axis=1)

        self.set_a(quant)
        self.scatter.set_offsets(coords[::self.po.a])
        self.scatter.set_array(quant[::self.po.a])

        if highlight_coords is not None:
            self.scatter.set_alpha(self.po.alpha_excluded)
            self.highlighted_scatter.set_alpha(self.po.alpha_included)
            highlight_coords = np.delete(highlight_coords, zc, axis=1)
            self.highlighted_scatter.set_offsets(highlight_coords[::self.po.a])
            self.highlighted_scatter.set_array(highlight_quant[::self.po.a])
            self.cb.set_alpha(1)
            self.cb.draw_all()

        self.set_title()
        self.set_scmp_lims()

    def update_plot_from_plot_options(self, plot_options):

        super(Scatter2D, self).update_plot_from_plot_options()

        self.po = plot_options
        tc, tq, hc, hq = self.compute_plot_content()
        self.update_plot(tc, tq, hc, hq)

    def compute_plot_content(self):

        super(Scatter2D, self).compute_plot_content()

        snap = self.po.ar.snapshots[self.po.t]
        coords = snap.get_from_h5(n.COORDINATES)
        quant = snap.get_from_h5(self.po.quantity)

        tc, tq = self.apply_location_cutoffs(coords, quant)
        tc, tq = self.apply_quantity_cutoffs(tc, tq)
        tc, tq, hc, hq = self.apply_highlight_cutoffs(tc, tq)

        # Take the magnitude of vector quantities
        if tq.ndim > 1:
            tq = np.apply_along_axis(np.linalg.norm, 1, tq)
            if hq is not None:
                hq = np.apply_along_axis(np.linalg.norm, 1, hq)

        return tc, tq, hc, hq
