from arepo_vis import make_pcolor, make_radial
from h5_file import ArepoSnapshot
from pyhelm_eos import loadhelm_eos
import matplotlib.pyplot as plt
from run import ArepoRun
from const import msol
from names import n
import numpy as np
import create_ics
import pprint
import ic


def test_make_radial():
    boxsize = 1e10
    a = [boxsize / 2, boxsize / 2, boxsize / 2]
    b = [boxsize / 1, boxsize / 2, boxsize / 2]
    cyl_rad = 0.01 * boxsize
    nshells = 200

    ar = ArepoRun.from_directory("/home/pierre/Desktop/AREPO/snapshots/"
                                 "97863140b478c1319cec0fd2f29258d6d36b5927/output_wd0_35He")
    s = ar.snapshots[1]

    quant = s.get_from_h5(n.DENSITY)
    coords = s.get_from_h5(n.COORDINATES)

    data = make_radial(coords.astype('float64'), quant.astype('float64'),
                       a, b,
                       cyl_rad,
                       nshells)
    fig, ax = plt.subplots()
    ax.plot(data[1, :], data[0, :], 'b')
    plt.show()


def test_make_pcolor():
    boxsize = 1e10
    resolutions = [1024, 1024]
    centers = [boxsize / 2, boxsize / 2, boxsize / 2]
    boxsizes = [1e10, 1e10]
    axes = [0, 1]

    ar = ArepoRun.from_directory("/home/pierre/Desktop/AREPO/snapshots/"
                                 "97863140b478c1319cec0fd2f29258d6d36b5927/output_wd0_35He")
    s = ar.snapshots[1]

    quant = s.get_from_h5(n.DENSITY)
    coords = s.get_from_h5(n.COORDINATES)

    data = make_pcolor(coords.astype('float64'), quant.astype('float64'),
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
    xnuc = [0.0, 0.5, 0.5, 0.0, 0.0]

    polytrope = ic.create_polytrope(eos, 3.0, 5e6, xnuc, pres_c=0.0, temp_c=5e5, dr=1e6)

    print(polytrope)


def test_create_wd():
    helm_file = "/home/pierre/PycharmProjects/arepo_helper/arepo_helper/data/eostable/helm_table.dat"
    species_file = "/home/pierre/PycharmProjects/arepo_helper/arepo_helper/data/eostable/species05.txt"
    eos = loadhelm_eos(helm_file, species_file, True)
    rho_c = 5e6
    xnuc = [0.0, 0.5, 0.5, 0.0, 0.0]
    wd = ic.create_wd(eos, rho_c, temp_c=5e5, xnuc_py=xnuc, tolerance=1e-6)

    print(wd)


def test_create_wd_wdec():
    # Generate 1d profile
    wdec_dir = "/home/pierre/wdec/"  # TODO: trailing foward slash is required!
    wd = ic.create_wd_wdec(wdec_dir, 5)

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

    e_calculated, dedt, p_calculated, csnd = eos.tgiven(rho_c, xnuc, t_c)
    print(f"{e_calculated=} {dedt=} {p_calculated=} {csnd=}")

    data_dict = eos.tgivenfull(rho_c, xnuc, t_c)
    pp.pprint(data_dict)

    data_dict = eos.ptgivenfull(pres_c, xnuc, t_c, rho_c)
    pp.pprint(data_dict)


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


def test_convert_to_healpix():
    # Generate 1d profile
    helm_file = "/home/pierre/PycharmProjects/arepo_helper/arepo_helper/data/eostable/helm_table.dat"
    species_file = "/home/pierre/PycharmProjects/arepo_helper/arepo_helper/data/eostable/species05.txt"
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
    helm_file = "/home/pierre/PycharmProjects/arepo_helper/arepo_helper/data/eostable/helm_table.dat"
    species_file = "/home/pierre/PycharmProjects/arepo_helper/arepo_helper/data/eostable/species05.txt"
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
                                                  boxres=32, grid_pres=2e6, grid_density=1e-4, grid_xnuc=xnuc)

    print(complete_dict)


def test_rho_c_from_mtot():
    mtot = 0.35  # msol
    temp_c = 5e5
    xnuc = [1.0, 0.0, 0.0, 0.0, 0.0]

    helm_file = "/home/pierre/PycharmProjects/arepo_helper/arepo_helper/data/eostable/helm_table.dat"
    species_file = "/home/pierre/PycharmProjects/arepo_helper/arepo_helper/data/eostable/species05.txt"
    eos = loadhelm_eos(helm_file, species_file, True)

    rho_c = ic.rho_c_from_mtot(mtot, temp_c, eos, xnuc)
    print(f"{rho_c=}")


def test_mtot_from_rho_c():
    rho_c = 1e4
    temp_c = 5e5
    xnuc = [1.0, 0.0, 0.0, 0.0, 0.0]

    helm_file = "/home/pierre/PycharmProjects/arepo_helper/arepo_helper/data/eostable/helm_table.dat"
    species_file = "/home/pierre/PycharmProjects/arepo_helper/arepo_helper/data/eostable/species05.txt"
    eos = loadhelm_eos(helm_file, species_file, True)

    mtot = ic.mtot_from_rho_c(rho_c, temp_c, eos, xnuc)
    print(f"{mtot}")


def test_create_particles_fill_grid():
    # Generate 1d profile
    helm_file = "/home/pierre/PycharmProjects/arepo_helper/arepo_helper/data/eostable/helm_table.dat"
    species_file = "/home/pierre/PycharmProjects/arepo_helper/arepo_helper/data/eostable/species05.txt"
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

    ugh = create_ics.create_particles_fill_grid(pos.astype('float64'), boxsize, 32)
    print(ugh)


def test_all_libs_functions():
    # PYHELM_EOS
    test_helm_eos()

    # IC
    test_create_wd_wdec()
    test_make_polytrope()
    test_create_wd()
    test_rho_c_from_mtot()
    test_mtot_from_rho_c()

    # CREATE_ICS
    test_convert_to_healpix()
    test_create_particles_fill_grid()
    test_add_grid_particles()

    # VISUALISE
    test_make_radial()
    test_make_pcolor()


def main():
    test_all_libs_functions()


if __name__ == '__main__':
    main()
