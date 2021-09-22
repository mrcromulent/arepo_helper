from arepo_helper.plot_manager import quick_pcolor, quick_radial, make_radial_time_series, quick_nuclear_pcolors
from arepo_helper.plot_manager import quick_nuclear_compositions, compute_plot_options_array
from arepo_helper.abstract_plot import mapping, get_plotter_func, get_plot_options_func
from arepo_helper.abstract_plot import PColorPlot, PColorPlotOptions
from arepo_helper.abstract_plot import RadialPlotOptions, RadialPlot
from arepo_helper.abstract_plot import ScatterPlotOptions, Scatter2D
from arepo_helper.group_animation import GroupAnimation
from arepo_helper.species import ArepoSpeciesList
from arepo_helper.analysis import ArepoAnalyser
from arepo_vis import make_pcolor, make_radial
from arepo_helper.utilities import SuppressStdout
from arepo_helper.run import ArepoRun
from pyhelm_eos import loadhelm_eos
from arepo_helper.const import msol
from arepo_helper.names import n
from definitions import DATA_DIR
import matplotlib.pyplot as plt
import numpy as np
import create_ics
import pyopal_eos
import pprint
import pytest
import pysph
import pyeos
import ic
import os


opal_file = os.path.join(DATA_DIR, "eostable", "EOS5_data")
helm_file = os.path.join(DATA_DIR, "eostable", "helm_table_541_201.dat")
species_file = os.path.join(DATA_DIR, "eostable", "species05.txt")
wdec_dir = os.path.join(DATA_DIR, "wdec")
arepo_dir = os.path.join(DATA_DIR, "snapshots")


class TestArepoHelper:
    pass


def test_utilities():
    from arepo_helper.utilities import Coordinates, get_cmap, SuppressStdout, sci, \
        convert_from_ruediger_dict, convert_to_ruediger_dict

    rue_dict = convert_to_ruediger_dict({n.DENSITY: None})
    convert_from_ruediger_dict(rue_dict)
    Coordinates.coordinates(["x", "y"])
    get_cmap(n.DENSITY, [1, 1e7])
    sci(100)
    with SuppressStdout():
        print(5)


def test_plot_manager():
    from arepo_helper.plot_manager import plot_particle_lifetimes
    ar = ArepoRun.from_directory(arepo_dir)
    asl = ArepoSpeciesList(species_file)
    s = ar.snapshots[0]

    quick_pcolor(s, n.DENSITY)
    quick_radial(s, n.DENSITY)
    quick_nuclear_compositions(s, asl)
    quick_nuclear_pcolors(s, asl)
    make_radial_time_series(ar, n.DENSITY)
    compute_plot_options_array(ar, [n.DENSITY, n.PRESSURE], [["x", "y"]],
                               plot_type=PColorPlotOptions.plot_type)
    plot_particle_lifetimes(ar, np.array([1, 2]))


def test_run():
    ar = ArepoRun.from_directory(arepo_dir)
    ar.get_created_destroyed_particles(1)
    ar.get_particle_lifetimes(np.array([1, 2]), [n.DENSITY])


def test_h5_file():
    from arepo_helper.h5_file import ArepoSnapshot
    asl = ArepoSpeciesList(species_file)
    ar = ArepoRun.from_directory(arepo_dir)
    center = np.array([5e9, 5e9, 5e9])
    s = ar.snapshots[0]
    s.get_from_h5(n.DENSITY)
    s.coords()
    s.mass()
    s.rho()
    s.vel()
    s.u()
    s.xnuc()
    s.pressure()
    s.r(center)
    s.min(n.DENSITY)
    s.max(n.DENSITY)
    s.get_field_names()
    s.mean_a(asl)
    s.center_of_mass()
    s.coords_center_at(center)
    s.angular_momentum(["x", "y"])
    s.num_particles()
    s.get_values_from_pids(n.DENSITY, np.array([1]))
    s.print_particle_info(np.array([1, 2]))
    s.get_ids_of_largest(n.DENSITY)
    s.get_ids_of_smallest(n.DENSITY)


def test_species():
    asl = ArepoSpeciesList(species_file)
    asl.index_of("He4")
    xnuc = np.array([1.0, 0.0, 0.0, 0.0, 0.0])
    asl.azbar(xnuc)
    asl.estimate_temp_from_e(1e16, xnuc)
    asl.estimate_e_from_temp(1e7, xnuc)


# @pytest.mark.skip(reason="Requires generalisation")
def test_simulation(tmpdir):
    from arepo_helper.sim_config import J, Paths
    from arepo_helper.h5_file import ArepoICs
    from arepo_helper.simulation import ArepoSimulation
    from create_ics import add_grid_particles
    import create_ics
    import ic

    # Top level project descriptor
    project_name = "default"

    simulation_out_dir = tmpdir

    # Other specifications
    arepo_version = "dissipative"
    xnuc = [1.0, 0.0, 0.0, 0.0, 0.0]
    wd_total_mass = 0.35
    pmass = 1e-6 * msol
    boxsize = 1e10
    time_max = 40
    temp_c = 5e5

    # Construct 1d profile
    asl = ArepoSpeciesList(species_file)
    eos = loadhelm_eos(helm_file, species_file, True)
    rho_c = ic.rho_c_from_mtot(wd_total_mass, temp_c, eos, xnuc)
    wd = ic.create_wd(eos, rho_c, temp_c, xnuc_py=xnuc)

    # Convert to healpix-distributed 3D particles
    healpix_return = create_ics.convert_to_healpix(wd, boxsize,
                                                   centers=[boxsize / 2, boxsize / 2, boxsize / 2],
                                                   randomizeshells=False,
                                                   randomizeradii=False,
                                                   pmass=pmass)

    # Add grid particles
    data = add_grid_particles(healpix_return, boxsize)

    # Meshrelax density
    data[n.MASSES] = data[n.DENSITY]

    # Jobscript
    js_explicit_options = {
        J.NAME: project_name,
        J.PROJECT_CODE: "y08",
        J.QUEUE_TYPE: "normal",
        J.WALLTIME: "23:59:00",
        J.EMAIL: "uri.burmester@anu.edu.au",
        J.MEMORY: "512gb",
        J.N_CPUS: "240",
        J.DIRECTORY: "wd",
    }

    # Config
    config_explicit_options = {

        # MHD
        "MHD": True,
        "MHD_SEEDFIELD": True,
        "RIEMANN_HLLD": True,

        # Refinement and derefinement
        "REFINEMENT_SPLIT_CELLS": True,
        "REFINEMENT_MERGE_CELLS": True,
        "REFINEMENT_VOLUME_LIMIT": True,

        # Mesh motion and regularization
        "REGULARIZE_MESH_CM_DRIFT": True,
        "REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED": True,
        "REGULARIZE_MESH_FACE_ANGLE": True,

        # Time integration options
        "TREE_BASED_TIMESTEPS": True,

        # Add or relax mesh
        "MESHRELAX_DENSITY_IN_INPUT": True,

        # Gravity treatment
        "GRAVITY_NOT_PERIODIC": True,

        # Gravity softening
        "ADAPTIVE_HYDRO_SOFTENING": True,
        "NSOFTTYPES_HYDRO": 64,

        # I/O
        "OUTPUT_PRESSURE": True,
        "EVALPOTENTIAL": True,
        "OUTPUTPOTENTIAL": True,
        "INPUT_IN_DOUBLEPRECISION": True,
        "OUTPUT_IN_DOUBLEPRECISION": True,

        # Degenerate Equation of State
        "EOS_NSPECIES": asl.num_species,

    }

    # Param
    param_explicit_options = {
        # Initial conditions
        "InitCondFile": f"{Paths.INPUT}/{project_name}",
        "MHDSeedDir": 0,
        "MHDSeedValue": 0,

        # Output file names and formats
        "OutputDir": f"{Paths.OUTPUT}",
        "EosTable": f"{Paths.INPUT}/{os.path.basename(helm_file)}",
        "EosSpecies": f"{Paths.INPUT}/{os.path.basename(species_file)}",
        "OutputListOn": 0,

        # Output frequency
        "TimeBetSnapshot": 1.0,
        "TimeBetStatistics": 0.1,
        "TimeOfFirstSnapshot": 0.0,

        # Simulated time span and spatial extent
        "BoxSize": boxsize,
        "PeriodicBoundariesOn": 0,
        "TimeBegin": 0.0,
        "TimeMax": time_max,

        # Cosmological parameters
        "ComovingIntegrationOn": 0,

        # Moving mesh
        "MaxVolume": 1e30,
        "MaxVolumeDiff": 10.0,
        "MinVolume": 0.0,
        "CellMaxAngleFactor": 2.25,

        # Refinement and derefinement
        "ReferenceGasPartMass": pmass,
        "TargetGasMassFactor": 1,
        "RefinementCriterion": 1,
        "DerefinementCriterion": 1,

        # Gravitational Softening
        "MinimumComovingHydroSoftening": 1e+06,
        "AdaptiveHydroSofteningSpacing": 1.2,

        # Cooling and star formation
        "CoolingOn": 0,
    }

    simulation = ArepoSimulation(
        project_name,
        simulation_out_dir,
        js_explicit_options,
        config_explicit_options,
        param_explicit_options,
        version=arepo_version)

    ic_name = os.path.join(tmpdir, project_name + ".hdf5")
    ic_file = ArepoICs(ic_name)
    ic_file.write_ics(data, simulation.config, simulation.param)

    simulation.validate_input_file(ic_name, helm_file, species_file)

    simulation.copy_files_to_input([
        species_file,
        helm_file,
        ic_name,
        __file__
    ])


def test_astro_utils():
    import arepo_helper.astro_utils as utils
    utils.estimate_collinear_lagrange_point(0.75, 0.35, 1)
    utils.test_p(1e7, 1e22)
    utils.p_wrt_e(1e7, 1e16)
    utils.e_wrt_p(1e7, 1e22)


def test_wdec(tmpdir):
    from arepo_helper.wdec_results import WDECResults

    wd = WDECResults(wdec_dir)
    wd.plot_results()
    wd.plot_abundances()
    wd.estimate_sound_crossing_time(helm_file, species_file)
    wd.write_temp_profile(os.path.join(tmpdir, "test.txt"))
    wd.get_xnuc_at_idx(0, None)
    wd.check_values(wd.gp, wd.ip)


def test_arepo_files():
    from arepo_helper.arepo_files import ArepoEnergyFile
    enf = ArepoEnergyFile(arepo_dir)
    enf.plot()


def test_pyeos():
    print(dir(pyeos))


def test_pyopal_eos():

    with SuppressStdout():
        eos = pyopal_eos.loadopal_eos(opal_file)
        print(dir(eos))
        print(eos.getrholim())
        print(eos.gettlim())
        print(eos.gettvalues())
        print(eos.getxlim())
        x_frac = 0.75
        rho = 6e7
        temp = 5e5
        print(eos.tgiven(x_frac, rho, temp))


@pytest.mark.skip(reason="Too long")
def test_pysph():

    ar = ArepoRun.from_directory(arepo_dir)
    s = ar.snapshots[0]
    coords = s.coords()
    masses = s.mass()

    with SuppressStdout():
        tree = pysph.makeTree(coords)

        #
        coord = np.array([5e9, 5e9, 5e9])
        hsml, density, weighted_neighbours = tree.calcHsml(coord, coords, masses, 4)
        print(hsml, density, weighted_neighbours)

        #
        dens = tree.calcDensity(coord, hsml, coords, masses)
        print(dens)

        #
        search_coords = np.array([[1e10, 1e10, 1e10], [5e9, 5e9, 5e9]], ndmin=2)
        numthreads = 1
        neighbours = tree.getNearestNeighbours(search_coords, coords, numthreads)
        print(neighbours)


def test_make_radial_time_series(tmpdir):
    ar = ArepoRun.from_directory(arepo_dir)

    bs  = 1e10
    a   = [bs / 2, bs / 2, bs / 2]
    b   = [bs / 1, bs / 2, bs / 2]
    with SuppressStdout():
        make_radial_time_series(ar, n.DENSITY, a=a, b=b)


def test_quick_nuclear_pcolors():
    asl = ArepoSpeciesList(species_file)
    ar = ArepoRun.from_directory(arepo_dir)
    s = ar.snapshots[1]
    quick_nuclear_pcolors(s, asl, {"inner_boxsize": 1e10})


def test_quick_nuclear_compositions():
    ar = ArepoRun.from_directory(arepo_dir)

    s = ar.snapshots[1]
    asl = ArepoSpeciesList(species_file)
    quick_nuclear_compositions(s, asl)


def test_quick_radial(tmpdir):
    ar = ArepoRun.from_directory(arepo_dir)

    s = ar.snapshots[1]
    plot = quick_radial(s, n.DENSITY)
    plot.save(os.path.join(tmpdir, "test.png"))


def test_quick_pcolor(tmpdir):
    ar = ArepoRun.from_directory(arepo_dir)

    s = ar.snapshots[1]
    plot = quick_pcolor(s, n.DENSITY)
    plot.save(os.path.join(tmpdir, "test.png"))


@pytest.mark.skip(reason="Too long")
def test_group_animation(tmpdir):

    ar = ArepoRun.from_directory(arepo_dir)

    for plot_type in mapping:
        poa = compute_plot_options_array(ar, [n.DENSITY, n.PRESSURE], [["x", "y"]],
                                         plot_type=plot_type)
        group_animation = GroupAnimation(poa, fps=1)
        group_animation.animate()
        group_animation.save(os.path.join(tmpdir, "test.mp4"))


@pytest.mark.skip(reason="Too long")
def test_animation(tmpdir):

    ar = ArepoRun.from_directory(arepo_dir)

    for plot_type in mapping:
        poa = compute_plot_options_array(ar, [n.DENSITY, n.PRESSURE], [["x", "y"]],
                                         plot_type=plot_type,
                                         explicit_options={"inner_boxsize": 5e9})
        animation = GroupAnimation(poa, fps=1)
        animation.animate()
        animation.save(os.path.join(tmpdir, "test.mp4"))


def test_frames(tmpdir):

    ar = ArepoRun.from_directory(arepo_dir)
    aa = ArepoAnalyser({"inner_boxsize": 5e9,
                        "t_idx": 0,
                        "quantity": n.DENSITY,
                        "orientation": ["x", "y"]
                        })

    for plot_type in mapping:
        plotter = get_plotter_func(plot_type)
        po_func = get_plot_options_func(plot_type)
        po = po_func(ar, aa)
        plot = plotter(po)
        plot.save(os.path.join(tmpdir, "test.png"))


def test_make_scatter_plot(tmpdir):

    ar = ArepoRun.from_directory(arepo_dir)
    aa_sp = ArepoAnalyser({"inner_boxsize": 1e10,
                           "cbar_lims": [1e-1, 1e7],
                           "log_cmap": True,
                           "quantity": n.DENSITY,
                           "orientation": ["x", "y"],
                           "t_idx": 1})

    spo = ScatterPlotOptions(ar, aa_sp)
    sp = Scatter2D(spo)
    sp.save(os.path.join(tmpdir, "test.png"))


def test_make_pcolor_plot(tmpdir):

    ar = ArepoRun.from_directory(arepo_dir)
    aa_pcp = ArepoAnalyser({"cbar_lims": [0.0, 1.0],
                            "log_cmap": False,
                            "quantity": n.NUCLEARCOMPOSITION,
                            "orientation": ["x", "y"],
                            "xlim": [0, 1e10],
                            "t_idx": 1,
                            "select_column": 0})

    pcpo = PColorPlotOptions(ar, aa_pcp)
    pcp = PColorPlot(pcpo)
    pcp.save(os.path.join(tmpdir, "test.png"))


def test_make_radial_plot(tmpdir):

    ar = ArepoRun.from_directory(arepo_dir)

    aa_rp = ArepoAnalyser({"inner_boxsize": 1e10,
                           "quantity": n.DENSITY,
                           "orientation": ["x", "y"],
                           "logscale": False,
                           "t_idx": 1})

    rpo = RadialPlotOptions(ar, aa_rp)
    rp = RadialPlot(rpo)
    rp.save(os.path.join(tmpdir, "test.png"))


def test_make_radial():
    boxsize = 1e10
    a = [boxsize / 2, boxsize / 2, boxsize / 2]
    b = [boxsize / 1, boxsize / 2, boxsize / 2]
    cyl_rad = 0.01 * boxsize
    nshells = 200

    ar = ArepoRun.from_directory(arepo_dir)
    s = ar.snapshots[1]

    quant = s.get_from_h5(n.DENSITY)
    coords = s.get_from_h5(n.COORDINATES)

    data = make_radial(coords.astype("float64"), quant.astype("float64"),
                       a, b,
                       cyl_rad,
                       nshells)
    fig, ax = plt.subplots()
    ax.plot(data[1, :], data[0, :], "b")
    plt.show()


def test_make_pcolor():
    boxsize = 1e10
    resolutions = [1024, 1024]
    centers = [boxsize / 2, boxsize / 2, boxsize / 2]
    boxsizes = [1e10, 1e10]
    axes = [0, 1]

    ar = ArepoRun.from_directory(arepo_dir)
    s = ar.snapshots[1]

    quant = s.get_from_h5(n.DENSITY)
    coords = s.get_from_h5(n.COORDINATES)

    data = make_pcolor(coords.astype("float64"), quant.astype("float64"),
                       axes,
                       boxsizes,
                       resolutions,
                       centers,
                       include_neighbours_in_output=1,
                       numthreads=1)

    fig, ax = plt.subplots()
    x = np.arange(resolutions[0] + 1, dtype="float64") / resolutions[0] * boxsizes[0] - 0.5 * boxsizes[0] + centers[0]
    y = np.arange(resolutions[1] + 1, dtype="float64") / resolutions[1] * boxsizes[1] - 0.5 * boxsizes[1] + centers[1]
    im = ax.pcolormesh(x, y, np.transpose(data["grid"]), shading="flat")
    plt.colorbar(im, ax=ax)
    plt.show()


def test_make_polytrope():

    eos = loadhelm_eos(helm_file, species_file, True)
    xnuc = [0.0, 0.5, 0.5, 0.0, 0.0]

    polytrope = ic.create_polytrope(eos, 3.0, 5e6, xnuc, pres_c=0.0, temp_c=5e5, dr=1e6)

    print(polytrope)


def test_create_wd():
    eos = loadhelm_eos(helm_file, species_file, True)
    rho_c = 5e6
    xnuc = [0.0, 0.5, 0.5, 0.0, 0.0]
    wd = ic.create_wd(eos, rho_c, temp_c=5e5, xnuc_py=xnuc, tolerance=1e-6)

    print(wd)


def test_create_wd_wdec():
    # Generate 1d profile
    wd = ic.create_wd_wdec(wdec_dir + "/", 5)  # TODO: trailing foward slash is required!

    # Convert to healpix-distributed 3D particles
    boxsize = 1e10
    centers = [boxsize / 2, boxsize / 2, boxsize / 2]
    randomizeshells = False
    randomizeradii = False
    pmass = 1e-6 * msol
    healpix_return = create_ics.convert_to_healpix(wd, boxsize,
                                                   centers=centers,
                                                   randomizeshells=randomizeshells,
                                                   randomizeradii=randomizeradii,
                                                   pmass=pmass)
    print(healpix_return)


def test_helm_eos():
    pp = pprint.PrettyPrinter(indent=4)
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

    e_calculated, dedt, p_calculated, csnd = eos.tgiven(rho_c, xnuc, t_c)
    print(f"{e_calculated=} {dedt=} {p_calculated=} {csnd=}")

    data_dict = eos.tgivenfull(rho_c, xnuc, t_c)
    pp.pprint(data_dict)

    data_dict = eos.ptgivenfull(pres_c, xnuc, t_c, rho_c)
    pp.pprint(data_dict)


def test_create_animation(tmpdir):

    ar = ArepoRun.from_directory(arepo_dir)

    poa = compute_plot_options_array(ar, [n.DENSITY, n.PRESSURE], [["x", "y"]],
                                     plot_type=PColorPlotOptions.plot_type)
    group_animation = GroupAnimation(poa, fps=1)
    group_animation.animate()
    group_animation.save(os.path.join(tmpdir, "test.mp4"))


def test_create_pcolor(tmpdir):
    ar = ArepoRun.from_directory(arepo_dir)
    aa = ArepoAnalyser({"inner_boxsize": 5e9,
                        "quantity": n.DENSITY,
                        "orientation": ["x", "y"]
                        })
    po = PColorPlotOptions(ar, aa)
    pcolor_plot = PColorPlot(po)
    pcolor_plot.save(os.path.join(tmpdir, "test.png"))


def test_convert_to_healpix():
    # Generate 1d profile
    eos = loadhelm_eos(helm_file, species_file, True)
    rho_c = 5e6
    xnuc = [0.0, 0.5, 0.5, 0.0, 0.0]
    wd = ic.create_wd(eos, rho_c, temp_c=5e5, xnuc_py=xnuc, tolerance=1e-6)

    # Convert to healpix-distributed 3D particles
    boxsize = 1e10
    centers = [boxsize / 2, boxsize / 2, boxsize / 2]
    randomizeshells = False
    randomizeradii = False
    pmass = 1e-6 * msol
    healpix_return = create_ics.convert_to_healpix(wd, boxsize,
                                                   centers=centers,
                                                   randomizeshells=randomizeshells,
                                                   randomizeradii=randomizeradii,
                                                   pmass=pmass)
    print(healpix_return)


def test_add_grid_particles():
    # Generate 1d profile
    eos = loadhelm_eos(helm_file, species_file, True)
    rho_c = 5e6
    xnuc = [1.0, 0.0, 0.0, 0.0, 0.0]
    wd = ic.create_wd(eos, rho_c, temp_c=5e5, xnuc_py=xnuc, tolerance=1e-6)

    # Convert to healpix-distributed 3D particles
    boxsize = 1e10
    centers = [boxsize / 2, boxsize / 2, boxsize / 2]
    randomizeshells = False
    randomizeradii = False
    pmass = 1e-6 * msol
    wd = create_ics.convert_to_healpix(wd, boxsize,
                                       centers=centers,
                                       randomizeshells=randomizeshells,
                                       randomizeradii=randomizeradii,
                                       pmass=pmass)

    complete_dict = create_ics.add_grid_particles(wd, boxsize,
                                                  boxres=32, grid_pres=4e6, grid_density=1e-4, grid_xnuc=xnuc)

    print(complete_dict)


def test_rho_c_from_mtot():
    mtot = 0.35  # msol
    temp_c = 5e5
    xnuc = [1.0, 0.0, 0.0, 0.0, 0.0]
    eos = loadhelm_eos(helm_file, species_file, True)

    rho_c = ic.rho_c_from_mtot(mtot, temp_c, eos, xnuc)
    print(f"{rho_c=}")


def test_mtot_from_rho_c():
    rho_c = 1e4
    temp_c = 5e5
    xnuc = [1.0, 0.0, 0.0, 0.0, 0.0]
    eos = loadhelm_eos(helm_file, species_file, True)

    mtot = ic.mtot_from_rho_c(rho_c, temp_c, eos, xnuc)
    print(f"{mtot}")


def test_create_particles_fill_grid():
    # Generate 1d profile
    eos = loadhelm_eos(helm_file, species_file, True)
    rho_c = 5e6
    xnuc = [0.0, 0.5, 0.5, 0.0, 0.0]
    wd = ic.create_wd(eos, rho_c, temp_c=5e5, xnuc_py=xnuc, tolerance=1e-6)

    # Convert to healpix-distributed 3D particles
    boxsize = 1e10
    centers = [boxsize / 2, boxsize / 2, boxsize / 2]
    randomizeshells = False
    randomizeradii = False
    pmass = 1e-6 * msol
    healpix_return = create_ics.convert_to_healpix(wd, boxsize,
                                                   centers=centers,
                                                   randomizeshells=randomizeshells,
                                                   randomizeradii=randomizeradii,
                                                   pmass=pmass)

    pos = healpix_return[n.COORDINATES]

    ugh = create_ics.create_particles_fill_grid(pos.astype("float64"), boxsize, 32)
    print(ugh)
