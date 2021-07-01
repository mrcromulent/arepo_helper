from snapshot import ArepoSnapshot
from pyhelm_eos import loadhelm_eos
import matplotlib.pyplot as plt
from const import msol
import arepo_radial
import arepo_pcolor
from names import n
import numpy as np
import create_ics
import pprint
import h5py
import ic


def test_temperature_variations():
    t = ["000", "001", "002", "003", "004", "005"]
    center = [5e9, 5e9, 5e9]

    for tval in t:
        filepath = f"/home/pierre/Desktop/3_healpix/output/snapshot_{tval}.hdf5"
        with h5py.File(filepath, 'r') as file:
            coordi = np.array(file["/PartType0/Coordinates"])
            tempi = np.array(file["/PartType0/Temperature"])
            densi = np.array(file["/PartType0/Density"])
            xnuci = np.array(file["/PartType0/NuclearComposition"])

            print("{:e}".format(np.average(coordi[:, 0])))
            print("{:e}".format(np.average(coordi[:, 1])))
            print("{:e}".format(np.average(coordi[:, 2])))

            i = (xnuci[::, 1] > 0) & (tempi > 1e3)
            selected_coords = coordi[i]
            selected_temps = tempi[i]
            r = np.sqrt(
                (selected_coords[:, 0] - center[0]) ** 2 +
                (selected_coords[:, 1] - center[1]) ** 2 +
                (selected_coords[:, 2] - center[2]) ** 2)
            inds = r.argsort()
            radius = r[inds]
            radial_temp = selected_temps[inds]

            fig, ax = plt.subplots()
            plt.plot(radius, radial_temp)
            ax.set_yscale('log')
            plt.show()


def test_temperature_with_density_cutoff():
    filepath = f"/home/pierre/Desktop/AREPO/AREPO_2020/arepo_new/snapshot_003.hdf5"
    center = [5e9, 5e9, 5e9]

    with h5py.File(filepath, 'r') as file:
        coordi = np.array(file["/PartType0/Coordinates"])
        densi = np.array(file["/PartType0/Density"])
        ui = np.array(file["/PartType0/InternalEnergy"])

        i = (densi > 1e-6)
        selected_coords = coordi[i]
        selected_temps = ui[i] * 1.3807e-16
        r = np.sqrt(
            (selected_coords[:, 0] - center[0]) ** 2 +
            (selected_coords[:, 1] - center[1]) ** 2 +
            (selected_coords[:, 2] - center[2]) ** 2)
        inds = r.argsort()
        radius = r[inds]
        radial_temp = selected_temps[inds]

        # print("{:e}".format(tempsum))
        fig, ax = plt.subplots()
        plt.plot(radius, radial_temp)
        ax.set_yscale('log')
        plt.show()


def test_make_radial():
    boxsize = 1e10
    a = np.array([boxsize / 2, boxsize / 2, boxsize / 2])
    b = np.array([boxsize, boxsize / 2, boxsize / 2])
    cyl_rad = 0.01 * boxsize
    nshells = 200

    s = ArepoSnapshot("/home/pierre/Desktop/AREPO/AREPO_2020/arepo_new/snapshot_000.hdf5")

    rho = s.get_from_h5(n.DENSITY)
    u = s.get_from_h5(n.INTERNALENERGY)
    pres = s.get_from_h5(n.PRESSURE)
    coords = s.get_from_h5(n.COORDINATES)

    data = arepo_radial.make_radial(coords.astype('float64'), pres.astype('float64'),
                                    a, b,
                                    cyl_rad,
                                    nshells)
    fig, ax = plt.subplots()
    ax.plot(data[1, :], data[0, :], 'b')
    plt.show()
    fig.savefig('test.png')


def test_make_pcolor():
    boxsize = 1e10
    resolutions = np.array([1024, 1024])
    centers = np.array([boxsize / 2, boxsize / 2, boxsize / 2])
    boxsizes = np.array([1e9, 1e9])
    axes = np.array([0, 1])

    s = ArepoSnapshot("/home/pierre/Desktop/test/snapshot_010.hdf5")

    rho = s.get_from_h5(n.DENSITY)
    u = s.get_from_h5(n.INTERNALENERGY)
    pres = s.get_from_h5(n.PRESSURE)
    coords = s.get_from_h5(n.COORDINATES)

    data = arepo_pcolor.make_pcolor(coords.astype('float64'), u.astype('float64'),
                                    axes,
                                    boxsizes,
                                    resolutions,
                                    centers,
                                    include_neighbours_in_output=1,
                                    numthreads=1)

    fig, ax = plt.subplots()
    x = np.arange(resolutions[0] + 1, dtype="float64") / resolutions[0] * boxsizes[0] - 0.5 * boxsizes[0] + centers[0]
    y = np.arange(resolutions[1] + 1, dtype="float64") / resolutions[1] * boxsizes[1] - 0.5 * boxsizes[1] + centers[1]
    im = ax.pcolormesh(x, y, np.transpose(data["grid"]), shading='flat')
    plt.colorbar(im, ax=ax)
    plt.show()


def test_make_polytrope():
    helm_file = "/home/pierre/PycharmProjects/arepo_helper/arepo_helper/data/eostable/helm_table.dat"
    species_file = "/home/pierre/PycharmProjects/arepo_helper/arepo_helper/data/eostable/species05.txt"
    eos = loadhelm_eos(helm_file, species_file, True)
    xnuc = np.array([0.0, 0.5, 0.5, 0.0, 0.0])

    polytrope = ic.create_polytrope(eos, 3.0, 5e6, xnuc, pres_c=0.0, temp_c=5e5, dr=1e6)

    print(polytrope)


def test_create_wd():
    helm_file = "/home/pierre/PycharmProjects/arepo_helper/arepo_helper/data/eostable/helm_table.dat"
    species_file = "/home/pierre/PycharmProjects/arepo_helper/arepo_helper/data/eostable/species05.txt"
    eos = loadhelm_eos(helm_file, species_file, True)
    rho_c = 5e6
    xnuc = np.array([0.0, 0.5, 0.5, 0.0, 0.0])
    wd = ic.create_wd(eos, rho_c, temp_c=5e5, xnuc_py=xnuc, tolerance=1e-6)

    print(wd)


def test_create_wd_wdec():
    # Generate 1d profile
    wdec_dir = "/home/pierre/wdec/"  # TODO: trailing foward slash is required!
    wd = ic.create_wd_wdec(wdec_dir)

    # Convert to healpix-distributed 3D particles
    nspecies = 5
    boxsize = 1e10
    centers = np.array([boxsize / 2, boxsize / 2, boxsize / 2])
    makebox = True
    randomizeshells = False
    randomizeradii = False
    pmass = 1e-6 * msol
    healpix_return = create_ics.convert_to_healpix(wd, nspecies, boxsize,
                                                   centers=centers,
                                                   makebox=makebox,
                                                   randomizeshells=randomizeshells,
                                                   randomizeradii=randomizeradii,
                                                   pmass=pmass)


def test_helm_eos():
    pp = pprint.PrettyPrinter(indent=4)

    helm_file = "/home/pierre/PycharmProjects/arepo_helper/arepo_helper/data/eostable/helm_table.dat"
    species_file = "/home/pierre/PycharmProjects/arepo_helper/arepo_helper/data/eostable/species05.txt"
    eos = loadhelm_eos(helm_file, species_file, True)
    xnuc = np.array([0.0, 0.5, 0.5, 0.0, 0.0])
    rho_c = 2e6
    pres_c = 1e23
    u_c = 8.5e16
    t_c = 6e8

    temp_calculated, u_calculated = eos.pgiven(rho_c, xnuc, pres_c)
    print(f"{temp_calculated=} {u_calculated=}")

    temp_calculated, p_calculated = eos.egiven(rho_c, xnuc, u_c)
    print(f"{temp_calculated=} {p_calculated=}")

    e_calculated, dedT, p_calculated, csnd = eos.tgiven(rho_c, xnuc, t_c)
    print(f"{e_calculated=} {dedT=} {p_calculated=} {csnd=}")

    data_dict = eos.tgivenfull(rho_c, xnuc, t_c)
    pp.pprint(data_dict)

    data_dict = eos.ptgivenfull(pres_c, xnuc, t_c, rho_c)
    pp.pprint(data_dict)


def test_check_mass():
    # s = ArepoSnapshot("/home/pierre/CLionProjects/arepo_helper_libs/cmake-build-debug/bin.dat.ic.hdf5")
    # s = ArepoSnapshot(
    # "/run/user/1000/gvfs/sftp:host=gadi.nci.org.au/g/data/y08/ub0692/simulations/singular_WDs/second_WD/output/snapshot_000.hdf5")
    for t in range(0, 10):
        fname = f"/home/pierre/Desktop/test/snapshot_00{t}.hdf5"
        s = ArepoSnapshot(fname)
        masses = s.get_from_h5(n.MASSES)
        density = s.get_from_h5(n.DENSITY)

        print(s)
        # print(f"{min(masses)=}")
        # print(f"{max(masses)=}")
        # print(f"{min(density)=}")
        # print(f"{max(density)=}")


def test_create_wd_python():
    from pyhelm_eos import loadhelm_eos
    from wd_utils import WDUtils

    equation_of_state = loadhelm_eos("./eostable/helm_table.dat", "./eostable/species05.txt", True)
    central_density = WDUtils.get_rho_c_from_mass(0.55, equation_of_state)
    print(f"Central Density = {'{:.5E}'.format(central_density)}")


def test_create_animation():
    from run import ArepoRun
    from analysis import ArepoAnalyser
    from plot_manager import PlotManager
    from group_animation import GroupAnimation

    # ar = ArepoRun.from_directory(
    # "/run/user/1000/gvfs/sftp:host=gadi.nci.org.au/scratch/y08/ub0692/test/dissipative_3/output")
    ar = ArepoRun.from_directory("/home/pierre/Desktop/output")
    aa = ArepoAnalyser(analysis_options={"inner_boxsize": 5e9})
    apm = PlotManager(ar, analyser=aa)

    poa = np.transpose(np.array((
        [apm.compute_plot_options(0, n.DENSITY, ["x", "y"], "Scatter", explicit_options={"log_cmap": True})]
    )))

    test = GroupAnimation([0, 20], poa, apm)
    test.animate()
    test.save()


def test_create_pcolor():
    from run import ArepoRun
    from analysis import ArepoAnalyser
    from plot_manager import PlotManager
    from pcolor_plot import PColorPlot
    # from matplotlib.cm import ScalarMappable
    # from matplotlib.colors import SymLogNorm
    from names import n

    ar = ArepoRun.from_directory("/home/pierre/Desktop/test")
    aa = ArepoAnalyser(analysis_options={"inner_boxsize": 1.3e9})
    apm = PlotManager(ar, analyser=aa)

    po = apm.compute_plot_options(8, n.INTERNALENERGY, ["x", "y"], "PColor", explicit_options={"log_cmap": True})

    pcolor_plot = PColorPlot(po)
    pcolor_plot.save()


def test_create_simulation():
    from simulation import ArepoSimulation
    from sim_config import J, Paths

    project_name = "dissipative"
    boxsize = 1e10

    js_explicit_options = {
        J.NAME: project_name,
        J.PROJECT_CODE: "y08",
        J.QUEUE_TYPE: "express",
        J.WALLTIME: "23:59:00",
        J.EMAIL: "uri.burmester@anu.edu.au",
        J.MEMORY: "512gb",
        J.N_CPUS: "240",
        J.DIRECTORY: "wd",
    }
    config_explicit_options = {

        "FORCE_EQUAL_TIMESTEPS": True,

        # MHD
        "MHD": True,
        "MHD_SEEDFIELD": True,
        "RIEMANN_HLLD": True,

        # EOS
        "EOS_NSPECIES": 5,

        #
        "REFINEMENT_SPLIT_CELLS": True,
        "REFINEMENT_MERGE_CELLS": True,
        "REFINEMENT_VOLUME_LIMIT": True,

        #
        "REGULARIZE_MESH_CM_DRIFT": True,
        "REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED": True,
        "REGULARIZE_MESH_FACE_ANGLE": True,

        # Relaxing
        "RELAXOBJECT": True,
        # "RELAXOBJECT_COOLING": True,

        #
        "GRAVITY_NOT_PERIODIC": True,

        # I/O
        "OUTPUT_PRESSURE": True,
        "INPUT_IN_DOUBLEPRECISION": True,
        "OUTPUT_IN_DOUBLEPRECISION": True
    }
    param_explicit_options = {
        # Initial conditions
        "InitCondFile": f"{Paths.INPUT}/bin.dat.ic",
        "MHDSeedDir": 0,
        "MHDSeedValue": 0,

        # Output file names and formats
        "OutputDir": f"{Paths.OUTPUT}",
        "EosTable": f"{Paths.INPUT}/helm_table.dat",
        "EosSpecies": f"{Paths.INPUT}/species05.txt",
        "OutputListOn": 0,

        # Output frequency
        "TimeBetSnapshot": 0.1,
        "TimeBetStatistics": 0.1,
        "TimeOfFirstSnapshot": 0.0,

        # Simulated time span and spatial extent
        "BoxSize": boxsize,
        "PeriodicBoundariesOn": 0,
        "TimeBegin": 0.0,
        "TimeMax": 1.0,
        "RelaxBaseFac": 0.01,
        # "RelaxTemperature": 5e5,

        # Cosmological parameters
        "ComovingIntegrationOn": 0,

        # Moving mesh
        "MaxVolume": boxsize ** 3,
        "MaxVolumeDiff": 10.0,
        "MinVolume": 0.0,
        "CellMaxAngleFactor": 2.25,
        # "CellShapingFactor": 1.0,

        # Refinement and derefinement
        "ReferenceGasPartMass": 2e27,
        "TargetGasMassFactor": 1,
        "RefinementCriterion": 1,
        "DerefinementCriterion": 1,

        # Cooling and star formation
        "CoolingOn": 0,
    }

    simulation = ArepoSimulation(
        project_name,
        "/home/pierre/",
        js_explicit_options,
        config_explicit_options,
        param_explicit_options,
        version=2019)

    simulation.copy_files_to_input([
        "/home/pierre/PycharmProjects/arepo_helper/arepo_helper/data/eostable/species05.txt",
        "/home/pierre/PycharmProjects/arepo_helper/arepo_helper/data/eostable/helm_table.dat",
        "/home/pierre/CLionProjects/arepo_helper_libs/cmake-build-debug/bin.dat.ic.hdf5"
    ])


def test_species():
    from species import ArepoSpeciesList
    a = ArepoSpeciesList("./snapshots/species05.txt")
    print(f"{a.num_species} found: {a.species_dict.keys()}")
    print(f"{a.species_dict['He4'].atomic_number}")


if __name__ == '__main__':
    # test_make_radial()
    # test_make_pcolor()
    # test_make_polytrope()
    # test_create_wd()
    # test_create_wd_wdec()
    # test_helm_eos()
    # test_check_mass()

    from run import ArepoRun
    from analysis import ArepoAnalyser
    from plot_manager import PlotManager
    from pcolor_plot import PColorPlot
    # from matplotlib.cm import ScalarMappable
    # from matplotlib.colors import SymLogNorm
    from utilities import suppress_stdout_stderr
    from utilities import Coordinates as C
    from names import n

    direc = "/home/pierre/Desktop/output/"
    inner_boxsize = 2e9
    tmax = 20

    ar = ArepoRun.from_directory(direc)
    aa = ArepoAnalyser(analysis_options={"inner_boxsize": inner_boxsize})
    apm = PlotManager(ar, analyser=aa)

    # for t in range(0, tmax):
    #
    #     po = apm.compute_plot_options(t, n.DENSITY, ["x", "y"], "PColor", explicit_options={"log_cmap": True})
    #     pl = PColorPlot(po)
    #     pl.save()

    # for t in range(0, tmax):
    #
    #     po = apm.compute_plot_options(t, n.DENSITY, ["x", "y"], "PColor", explicit_options={"log_cmap": True})
    #
    #     s = ArepoSnapshot(f"{direc}snapshot_{t:03}.hdf5")
    #     coords = s.get_from_h5(n.COORDINATES)
    #     quant = s.get_from_h5(n.DENSITY)
    #     res = 1204
    #     boxsize_x = po.xlim[1] - po.xlim[0]
    #     boxsize_y = po.ylim[1] - po.ylim[0]
    #     resolutions = np.array([res, res])
    #     boxsizes = np.array([boxsize_x, boxsize_y])
    #     centers = np.array([np.average(po.xlim),
    #                         np.average(po.ylim),
    #                         np.average(po.zlim)])
    #     xc, yc, _ = C.coordinates(po.orientation)
    #     axes = np.array([xc, yc])
    #
    #     with suppress_stdout_stderr():
    #         data = arepo_pcolor.make_pcolor(coords.astype('float64'), quant.astype('float64'),
    #                                         axes,
    #                                         boxsizes,
    #                                         resolutions,
    #                                         centers,
    #                                         include_neighbours_in_output=1,
    #                                         numthreads=1)
    #
    #     x = np.arange(res + 1, dtype="float64") / res * boxsize_x - 0.5 * boxsize_x + centers[0]
    #     y = np.arange(res + 1, dtype="float64") / res * boxsize_y - 0.5 * boxsize_y + centers[1]
    #
    #     fig, ax = plt.subplots()
    #     pcolor = ax.pcolormesh(x, y, np.transpose(data["grid"]), shading='flat')
    #     fig.colorbar(pcolor)
    #     filename = f"./{po.t}_{po.quantity}_{po.orientation}.png"
    #     fig.savefig(filename, dpi=300)

    radials = []
    labels = []
    for t in range(19,20):
        s = ArepoSnapshot(f"{direc}snapshot_{t:03}.hdf5")
        coords = s.get_from_h5(n.COORDINATES)
        quant = s.get_from_h5(n.NUCLEARCOMPOSITION)[:, 0]

        boxsize = 1e10
        a = np.array([boxsize / 2, boxsize / 2, boxsize / 2])
        b = np.array([boxsize, boxsize / 2, boxsize / 2])
        cyl_rad = 0.1 * boxsize
        nshells = 200
        radial = arepo_radial.make_radial(coords.astype('float64'),
                                          quant.astype('float64'),
                                          a, b,
                                          cyl_rad,
                                          nshells)

        radials.append(radial)
        labels.append(f"t = {round(s.get_from_h5(n.TIME))} s")

    fig, ax = plt.subplots()
    lines = []

    for i, radial in enumerate(radials):
        line = ax.plot(radial[1, :], radial[0, :], label=f"{labels[i]}")
        # ax.set_yscale('log')
        lines.append(line)

    plt.xlabel("Radial distance [cm]")
    plt.ylabel("xHe")
    plt.legend()
    plt.show()
    fig.savefig(f'radial_he_19.png')
