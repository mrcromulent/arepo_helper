from abstract_plot import PlotOptions, PColorPlot, PColorPlotOptions
from abstract_plot import RadialPlotOptions, RadialPlot, get_plot_options_func
from analysis import ArepoAnalyser
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from run import ArepoRun
from names import n
import numpy as np
import os


def quick_pcolor(s, quantity_name, inner_boxsize=None):
    ar = ArepoRun([s])
    aa = ArepoAnalyser({"inner_boxsize": inner_boxsize,
                        "quantity": quantity_name,
                        "orientation": ["x", "y"],
                        "t_idx": 0
                        })

    po = PColorPlotOptions(ar, aa)
    return PColorPlot(po)


def quick_radial(s, quantity_name, a=None, b=None, logscale=False, inner_boxsize=None):

    ar = ArepoRun([s])
    aa = ArepoAnalyser({"a": a,
                        "b": b,
                        "logscale": logscale,
                        "quantity": quantity_name,
                        "t_idx": 0,
                        "orientation": ["x", "y"],  # TODO: Refactor to exclude
                        "inner_boxsize": inner_boxsize,
                        })

    po = RadialPlotOptions(ar, aa)
    return RadialPlot(po)


def quick_nuclear_compositions(s, asl, a=None, b=None, logscale=False):
    ar = ArepoRun([s])
    lines   = []
    fig, ax = plt.subplots()

    for i, spec in enumerate(asl.species_dict):
        aa = ArepoAnalyser({"a": a,
                            "b": b,
                            "logscale": logscale,
                            "quantity": n.NUCLEARCOMPOSITION,
                            "select_column": i,
                            "t_idx": 0,
                            "orientation": ["x", "y"],  # TODO: Refactor to exclude
                            })

        po = RadialPlotOptions(ar, aa)
        rp = RadialPlot(po, figure=fig, ax=ax, do_not_plot=True)

        lines.append((rp.x, rp.y, spec))

    xnuc_colours = {
        "He4": "b",
        "C12": "k",
        "O16": "r",
        "Ne20": "g",
        "Mg24": "grey"
    }

    for i, radial in enumerate(lines):
        spec = radial[2]
        ax.plot(radial[0], radial[1], xnuc_colours[spec], label=f"{spec}")

    plt.ylim((0, 1.1))
    plt.xlabel("Radial distance [cm]")
    plt.ylabel(f"Nuclear Composition")
    plt.legend()
    plt.show()

    return fig


def quick_nuclear_pcolors(s, asl, inner_boxsize=None):

    ar = ArepoRun([s])
    for i, spec in enumerate(asl.species_dict):

        print(f"Creating PColor for species {spec}")

        aa = ArepoAnalyser({"inner_boxsize": inner_boxsize,
                            "quantity": n.NUCLEARCOMPOSITION,
                            "orientation": ["x", "y"],
                            "t_idx": 0,
                            "select_column": i,
                            "title": f"Location of {spec}",
                            "log_cmap": False,
                            })

        po = PColorPlotOptions(ar, aa)
        pcp = PColorPlot(po)
        pcp.save(f"{i}_{spec}.png")


def make_radial_time_series(ar, quantity_name, a=None, b=None, logscale=False):

    fig, ax = plt.subplots()
    nsnaps  = len(ar.snapshots)
    colors  = cm.viridis(np.linspace(0, 1, nsnaps))
    rp      = None

    for i, s in enumerate(ar.snapshots):
        aa = ArepoAnalyser({"a": a,
                            "b": b,
                            "logscale": logscale,
                            "quantity": quantity_name,
                            "t_idx": i,
                            "orientation": ["x", "y"],  # TODO: Refactor to exclude
                            "color": colors[i],
                            })

        po = RadialPlotOptions(ar, aa)
        rp = RadialPlot(po, figure=fig, ax=ax)
        t = s.get_from_h5(n.TIME)
        rp.line[0].set_label(f"t = {round(t, 2)}")

    #
    ax.legend()
    plt.title("Radial time series")
    rp.save()


def energy_balance(ar, export_filename=None):

    ie_list = []
    ke_list = []
    gpe_list = []
    t_list = []

    for i, s in enumerate(ar.snapshots):
        mass = s.get_from_h5(n.MASSES)
        velocity = s.get_from_h5(n.VELOCITIES)
        potential = s.get_from_h5(n.POTENTIAL)
        bfield = s.get_from_h5(n.MAGNETICFIELD)

        assert (np.all(bfield == 0.0))

        ie = mass * s.get_from_h5(n.INTERNALENERGY)
        ke = 0.5 * mass * np.linalg.norm(velocity, axis=1) ** 2
        gpe = mass * np.abs(potential)

        ie_list.append(ie.sum())
        ke_list.append(ke.sum())
        gpe_list.append(gpe.sum())
        t_list.append(s.get_from_h5(n.TIME))

    ke_list = np.array(ke_list)
    gpe_list = np.array(gpe_list)
    ie_list = np.array(ie_list)

    energy_sum = ie_list + ke_list - gpe_list

    fig, ax = plt.subplots()
    ax.plot(t_list, ie_list, label="InternalEnergy")
    ax.plot(t_list, ke_list, label="KineticEnergy")
    ax.plot(t_list, -gpe_list, label="Gravitational Potential Energy")
    ax.plot(t_list, energy_sum, label="Energy Sum")
    plt.xlabel("Simulation time [s]")
    plt.ylabel("Energy [erg]")
    plt.legend()

    max_val = max(energy_sum)
    min_val = min(energy_sum)
    max_pc_diff = abs((max_val - min_val) / min_val)
    print(f"Energy p/c diff: {max_pc_diff * 100}")

    if export_filename is not None:
        dir_name = os.path.dirname(ar.snapshots[0].filename)
        fn = os.path.join(dir_name, export_filename)
        fig.savefig(fn, dpi=300)
    else:
        plt.show()


def compute_plot_options_array(ar, qs, oris,
                               plot_type=PColorPlotOptions.plot_type,
                               explicit_options=None, ts=None):
    """
    :param ar: Arepo run
    :param qs: list of quantities
    :param oris: list of orientations
    :param plot_type: string of plot type
    :return: np.array of ArepoPlotOptions objects
    """

    if explicit_options is None:
        explicit_options = dict()

    if ts is None:
        ts = list(range(0, len(ar.snapshots)))

    poa     = np.empty((len(oris), len(qs), len(ts)), dtype=PlotOptions)
    po_func = get_plot_options_func(plot_type)

    for i, o in enumerate(oris):
        for j, q in enumerate(qs):
            for k, t in enumerate(ts):

                analysis_options = {"quantity": q,
                                    "orientation": o,
                                    "t_idx": t
                                    }
                analysis_options.update(explicit_options)
                aa = ArepoAnalyser(analysis_options)

                poa[i, j, k] = po_func(ar, aa)

    return poa
