from snapshot import ArepoSnapshot
from pyhelm_eos import loadhelm_eos
import matplotlib.pyplot as plt
import arepo_radial
import arepo_pcolor
from names import n
import numpy as np
import calcGrid
import create_ics
from const import msol
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


def test_calcRadialProfile_old():
    s = ArepoSnapshot("/home/pierre/Desktop/AREPO/AREPO_2020/arepo_new/snapshot_000.hdf5")
    c = np.array([5e9, 5e9, 5e9])

    u = s.get_from_h5(n.INTERNALENERGY)
    rho = s.get_from_h5(n.DENSITY)
    pres = s.get_from_h5(n.PRESSURE)
    coords = s.get_from_h5(n.COORDINATES)

    data = calcGrid.calcRadialProfile(coords.astype('float64'),
                                      u.astype('float64'),
                                      2,
                                      200,
                                      0,
                                      *c)

    fig, ax = plt.subplots()
    ax.plot(data[1, :], data[0, :], 'b')
    plt.show()


def test_calcASlice_old():
    res = 500
    c = np.array([5e9, 5e9, 5e9])
    ibsx = 1e10
    ibsy = 1e10

    s = ArepoSnapshot("/home/pierre/Desktop/AREPO/AREPO_2020/arepo_new/snapshot_000.hdf5")

    rho = s.get_from_h5(n.DENSITY)
    u = s.get_from_h5(n.INTERNALENERGY)
    pres = s.get_from_h5(n.PRESSURE)
    coords = s.get_from_h5(n.COORDINATES)

    data = calcGrid.calcASlice(coords.astype('float64'), u.astype('float64'),
                               res, res,
                               ibsx, ibsy,
                               *c,
                               0, 1,
                               proj=False,
                               boxz=0,
                               nz=1,
                               numthreads=1)

    fig, ax = plt.subplots()
    x = np.arange(res + 1, dtype="float64") / res * ibsx - 0.5 * ibsx + c[0]
    y = np.arange(res + 1, dtype="float64") / res * ibsy - 0.5 * ibsy + c[1]
    im = ax.pcolormesh(x, y, np.transpose(data["grid"]), shading='flat')
    plt.colorbar(im, ax=ax)
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


def test_make_pcolor():
    boxsize = 1e10
    resolutions = np.array([1024, 1024])
    centers = np.array([boxsize / 2, boxsize / 2, boxsize / 2])
    boxsizes = np.array([boxsize, boxsize])
    axes = np.array([0, 1])

    s = ArepoSnapshot("/home/pierre/Desktop/AREPO/AREPO_2020/arepo_new/snapshot_000.hdf5")

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

    helm_file       = "/home/pierre/PycharmProjects/arepo_helper/arepo_helper/data/eostable/helm_table.dat"
    species_file    = "/home/pierre/PycharmProjects/arepo_helper/arepo_helper/data/eostable/species05.txt"
    eos             = loadhelm_eos(helm_file, species_file, True)
    xnuc            = np.array([0.0, 0.5, 0.5, 0.0, 0.0])
    rho_c           = 2e6
    pres_c          = 1e23
    u_c             = 8.5e16
    t_c             = 6e8

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


if __name__ == '__main__':
    test_make_radial()
    test_make_pcolor()
    test_make_polytrope()
    test_create_wd()
    test_create_wd_wdec()
    test_helm_eos()
