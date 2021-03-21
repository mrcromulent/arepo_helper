class PlotOptions(object):

    plot_type   = None

    t = 0
    ar = None
    aa = None
    quantity = None
    orientation = None
    xlim = [None, None]
    ylim = [None, None]
    title = ""
    xlabel = ""
    ylabel = ""
    include_colorbar = True
    cbar_lims = [None, None]

    def __init__(self, *args, **kwargs):
        for dictionary in args:
            for key in dictionary:
                setattr(self, key, dictionary[key])
        for key in kwargs:
            setattr(self, key, kwargs[key])


class ScatterPlotOptions(PlotOptions):

    plot_type   = "Scatter"

    psize           = 0.01
    a               = 2
    alpha_excluded  = 0.05
    alpha_included  = 1
    cmap            = None
    scmp            = None

    def __init__(self, *args, **kwargs):
        super(ScatterPlotOptions, self).__init__(*args, **kwargs)


class RadialPlotOptions(PlotOptions):

    plot_type   = "Radial"

    logscale            = True
    nshells             = 200
    dr                  = 0
    color               = "k"
    include_colorbar    = False

    def __init__(self, *args, **kwargs):
        super(RadialPlotOptions, self).__init__(*args, **kwargs)


class PColorPlotOptions(PlotOptions):

    plot_type   = "PColor"

    resolution  = 1024
    projection  = False
    proj_fact   = 0.5
    boxz        = 0
    numthreads  = 1

    def __init__(self, *args, **kwargs):
        super(PColorPlotOptions, self).__init__(*args, **kwargs)


mapping = {ScatterPlotOptions.plot_type:   ScatterPlotOptions,
           PColorPlotOptions.plot_type:    PColorPlotOptions,
           RadialPlotOptions.plot_type:    RadialPlotOptions}
