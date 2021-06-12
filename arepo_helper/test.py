from snapshot import ArepoSnapshot
import matplotlib.pyplot as plt
import radial_pierre
import pcolor_pierre
from names import n
import numpy as np
import calcGrid
import h5py


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
    filepath = f"/home/pierre/Desktop/arepo_new/snapshot_003.hdf5"
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
    s = ArepoSnapshot("/home/pierre/Desktop/arepo_new/snapshot_000.hdf5")
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


def test_make_radial():
    res = 500
    c = np.array([5e9, 5e9, 5e9])
    ibsx = 1e10
    ibsy = 1e10

    s = ArepoSnapshot("/home/pierre/Desktop/arepo_new/snapshot_000.hdf5")

    rho = s.get_from_h5(n.DENSITY)
    u = s.get_from_h5(n.INTERNALENERGY)
    pres = s.get_from_h5(n.PRESSURE)
    coords = s.get_from_h5(n.COORDINATES)

    data = radial_pierre.make_radial(coords.astype('float64'), pres.astype('float64'), 200)
    fig, ax = plt.subplots()
    ax.plot(data[1, :], data[0, :], 'b')
    plt.show()


def test_calcASlice_old():
    res = 500
    c = np.array([5e9, 5e9, 5e9])
    ibsx = 1e10
    ibsy = 1e10

    s = ArepoSnapshot("/home/pierre/Desktop/arepo_new/snapshot_000.hdf5")

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


def test_make_pcolor():
    res = 500
    c = np.array([5e9, 5e9, 5e9])
    ibsx = 1e10
    ibsy = 1e10

    s = ArepoSnapshot("/home/pierre/Desktop/arepo_new/snapshot_000.hdf5")

    rho = s.get_from_h5(n.DENSITY)
    u = s.get_from_h5(n.INTERNALENERGY)
    pres = s.get_from_h5(n.PRESSURE)
    coords = s.get_from_h5(n.COORDINATES)

    data = pcolor_pierre.make_pcolor(coords.astype('float64'), u.astype('float64'),
                                     res, res,
                                     ibsx, ibsy, 0,
                                     *c,
                                     axis0=0, axis1=1,
                                     nz=1,
                                     include_neighbours_in_output=1,
                                     numthreads=1)

    fig, ax = plt.subplots()
    x = np.arange(res + 1, dtype="float64") / res * ibsx - 0.5 * ibsx + c[0]
    y = np.arange(res + 1, dtype="float64") / res * ibsy - 0.5 * ibsy + c[1]
    im = ax.pcolormesh(x, y, np.transpose(data["grid"]), shading='flat')
    plt.colorbar(im, ax=ax)
    plt.show()


if __name__ == '__main__':
    test_make_radial()
