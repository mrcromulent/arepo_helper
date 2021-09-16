from abstract_plot import PlotOptions, PColorPlot, PColorPlotOptions
from abstract_plot import RadialPlotOptions, RadialPlot, get_plot_options_func
from species import ArepoSpeciesList
from analysis import ArepoAnalyser
from h5_file import ArepoSnapshot
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from run import ArepoRun
from typing import Any
from names import n
import numpy as np
import os


def quick_pcolor(s: ArepoSnapshot,
                 quantity_name: str,
                 explicit_options: dict = None) -> PColorPlot:
    """Shortcut to generate a PColorPlot

    :param s: Snapshot
    :type s: ArepoSnapshot
    :param quantity_name: String name of quantity to be examined
    :type quantity_name: str
    :param explicit_options: Explicit options
    :type explicit_options: dict
    :return: PColorPlot
    :rtype: PColorPlot
    """

    if explicit_options is None:
        explicit_options = dict()

    analysis_options = {"quantity": quantity_name,
                        "orientation": ["x", "y"],
                        "t_idx": 0
                        }

    analysis_options.update(explicit_options)

    ar = ArepoRun([s])
    aa = ArepoAnalyser(analysis_options)

    po = PColorPlotOptions(ar, aa)
    return PColorPlot(po)


def quick_radial(s: ArepoSnapshot,
                 quantity_name: str,
                 a: list[float] = None,
                 b: list[float] = None,
                 explicit_options: dict = None) -> RadialPlot:
    """Shortcut to create RadialPlot

    :param s: Snapshot
    :type s: ArepoSnapshot
    :param quantity_name: Quantity name
    :type quantity_name: str
    :param a: Coordinates of start of radial cylinder
    :type a: list
    :param b: Coordinates of end of radial cylinder
    :type b: list
    :param explicit_options: Explicit options
    :type explicit_options: dict

    :return: RadialPlot
    :rtype: RadialPlot
    """

    if explicit_options is None:
        explicit_options = dict()

    analysis_options = {"a": a,
                        "b": b,
                        "logscale": False,
                        "quantity": quantity_name,
                        "t_idx": 0,
                        "orientation": ["x", "y"],  # TODO: Refactor to exclude
                        }
    analysis_options.update(explicit_options)

    ar = ArepoRun([s])
    aa = ArepoAnalyser(analysis_options)

    po = RadialPlotOptions(ar, aa)
    return RadialPlot(po)


def quick_nuclear_compositions(s: ArepoSnapshot,
                               asl: ArepoSpeciesList,
                               a: list[float] = None,
                               b: list[float] = None,
                               explicit_options: dict = None) -> Any:
    """Shortcut to create a radial plot showing location of different species.

    :param s: Snapshot
    :type s: ArepoSnapshot
    :param asl: Arepo species list object
    :type asl: ArepoSpeciesList
    :param a: Coordinates of start of radial cylinder
    :type a: list
    :param b: Coordinates of end of radial cylinder
    :type b: list
    :param explicit_options: Explicit options
    :type explicit_options: dict

    :return: Matplotlib figure
    :rtype: matplotlib.figure.Figure
    """

    if explicit_options is None:
        explicit_options = dict()

    ar = ArepoRun([s])
    lines = []
    fig, ax = plt.subplots()

    for i, spec in enumerate(asl.species_dict):
        analysis_options = {"a": a,
                            "b": b,
                            "logscale": False,
                            "quantity": n.NUCLEARCOMPOSITION,
                            "select_column": i,
                            "t_idx": 0,
                            "orientation": ["x", "y"],  # TODO: Refactor to exclude
                            }
        analysis_options.update(explicit_options)
        aa = ArepoAnalyser(analysis_options)

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


def quick_nuclear_pcolors(s: ArepoSnapshot,
                          asl: ArepoSpeciesList,
                          explicit_options: dict = None) -> None:
    """Shortcut to create a pcolor plots showing location of different species.

    :param s: Snapshot
    :type s: ArepoSnapshot
    :param asl: Arepo species list object
    :type asl: ArepoSpeciesList
    :param explicit_options: Explicit options
    :type explicit_options: dict

    """

    if explicit_options is None:
        explicit_options = dict()

    ar = ArepoRun([s])
    for i, spec in enumerate(asl.species_dict):
        print(f"Creating PColor for species {spec}")

        analysis_options = {"quantity": n.NUCLEARCOMPOSITION,
                            "orientation": ["x", "y"],
                            "t_idx": 0,
                            "select_column": i,
                            "title": f"Location of {spec}",
                            "log_cmap": False,
                            }

        analysis_options.update(explicit_options)
        aa = ArepoAnalyser(analysis_options)

        po = PColorPlotOptions(ar, aa)
        pcp = PColorPlot(po)
        pcp.save(f"{i}_{spec}.png")


def make_radial_time_series(ar: ArepoRun,
                            quantity_name: str,
                            a: list[float] = None,
                            b: list[float] = None,
                            explicit_options: dict = None) -> None:
    """Makes a plot showing radials across a number of snapshots

    :param ar: Arepo Run
    :type ar: ArepoRun
    :param quantity_name: Name of relevant quantity
    :type quantity_name: str
    :param a: Coordinates of start of radial cylinder
    :type a: list
    :param b: Coordinates of end of radial cylinder
    :type b: list
    :param explicit_options: Explicit options
    :type explicit_options: dict
    """

    if explicit_options is None:
        explicit_options = dict()

    fig, ax = plt.subplots()
    nsnaps = len(ar.snapshots)
    colors = cm.viridis(np.linspace(0, 1, nsnaps))
    rp = None

    for i, s in enumerate(ar.snapshots):

        analysis_options = {"a": a,
                            "b": b,
                            "logscale": False,
                            "quantity": quantity_name,
                            "t_idx": i,
                            "orientation": ["x", "y"],  # TODO: Refactor to exclude
                            "color": colors[i],
                            }
        analysis_options.update(explicit_options)
        aa = ArepoAnalyser(analysis_options)

        po = RadialPlotOptions(ar, aa)
        rp = RadialPlot(po, figure=fig, ax=ax)
        t = float(s.get_from_h5(n.TIME))
        rp.line[0].set_label(f"t = {round(t, 2)}")

    #
    ax.legend()
    plt.title(f"Radial time series: {quantity_name}")
    rp.save(f"Radial_Time_Series_n={nsnaps}_{quantity_name}.png")


def energy_balance(ar: ArepoRun,
                   export_filename: str = None) -> None:
    """Checks energy balance of all the energy types in an Arepo Run

    :param ar: Arepo Run
    :type ar: ArepoRun
    :param export_filename: Location to save plot to
    :type export_filename: str
    """

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


def compute_plot_options_array(ar: ArepoRun,
                               qs: list[str],
                               oris: list,
                               plot_type: str = PColorPlotOptions.plot_type,
                               explicit_options: dict = None,
                               ts: list[int] = None) -> np.ndarray:
    """Computes an array of plot options with (orientations, quantities, time) representing axes (0,1,2) respectively

    :param ar: Arepo run
    :type ar: ArepoRun
    :param qs: list of quantities
    :type qs: list
    :param oris: list of orientations
    :type oris: list
    :param plot_type: string of plot type
    :type plot_type: str
    :param explicit_options: Explicit options
    :type explicit_options: dict
    :param ts: Explicit specification of the relevant time indices
    :type ts: list

    :return: np.array of ArepoPlotOptions objects
    :rtype: np.ndarray
    """

    if explicit_options is None:
        explicit_options = dict()

    if ts is None:
        ts = list(range(0, len(ar.snapshots)))

    poa = np.empty((len(oris), len(qs), len(ts)), dtype=PlotOptions)
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


def plot_particle_lifetimes(ar: ArepoRun,
                            list_of_pids: list[int],
                            limit: int = 10) -> None:

    if len(list_of_pids) > limit:
        list_of_pids = list_of_pids[:limit]

    data = np.empty((len(list_of_pids), len(ar.snapshots), 3))

    for t_idx, s in enumerate(ar.snapshots):
        for i, p in enumerate(list_of_pids):
            if s.particle_exists(p):
                data[i, t_idx, :] = s.get_values_from_pids(n.COORDINATES, np.array([p]))
            else:
                data[i, t_idx, :] = np.nan

    fig, ax = plt.subplots()
    for i, p in enumerate(list_of_pids):
        ax.scatter(data[i, :, 0], data[i, :, 1],
                   c=data[i, :, 1],
                   s=1.5)
    for i, p in enumerate(list_of_pids):
        ax.annotate(f"{p}", (data[i, 0, 0], data[i, 0, 1]))

    fig.show()
