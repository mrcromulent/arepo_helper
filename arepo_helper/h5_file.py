"""
Tried to use hdfdict but it flat out didn't work. Did not dump files correctly nor did it preserve attributes of groups
"""

from utilities import Coordinates as C, suppress_stdout_stderr, part_fields
from arepo_vis import make_pcolor, make_radial
from names import ArepoHeader as h, n, path
import matplotlib.pyplot as plt
from pcolor_plot import PColorPlot
from radial_plot import RadialPlot
from plot_options import PColorPlotOptions, RadialPlotOptions
from matplotlib import colors
import numpy as np
import utilities
import h5py
import os


class ArepoH5File(object):

    def __init__(self, filename):
        self.data = dict()
        self.maxima = dict()
        self.minima = dict()
        self.filename = os.path.abspath(filename)

    def r(self, center):
        # TODO: Adapt for lower dimensionality

        coords = self.get_from_h5(n.COORDINATES)

        return np.sqrt(
            (coords[:, 0] - center[0]) ** 2 +
            (coords[:, 1] - center[1]) ** 2 +
            (coords[:, 2] - center[2]) ** 2)

    def center_of_mass(self):
        # TODO: Adapt for lower dimensionality

        coords = self.get_from_h5(n.COORDINATES)
        mass = self.get_from_h5(n.MASSES)

        rcm = [(coords[:, i] / mass.sum() * mass).sum() for i in [0, 1, 2]]
        return np.array(rcm)

    def mean_a(self, species_list):

        xnuc = self.get_from_h5(n.NUCLEARCOMPOSITION)

        mean_a = np.zeros(xnuc.shape[0])
        for i in range(species_list.num_species):
            a = list(species_list.species_dict.values())
            mean_a += xnuc[:, i] * a[i].mass_number
        return mean_a

    def coords_center_at(self, center):
        coords = self.get_from_h5(n.COORDINATES)
        return coords - np.array(center)[None, :]

    def angular_momentum(self):
        # TODO: this needs to be generalised for different orientations

        velocity = self.get_from_h5(n.VELOCITIES)
        mass = self.get_from_h5(n.MASSES)
        rcm = self.center_of_mass()
        coords = self.coords_center_at(rcm)

        angmom = (mass.astype(np.float64) *
                  (+ (coords.astype(np.float64)[:, 0] - rcm[0]) * velocity.astype(np.float64)[:, 1]
                   - (coords.astype(np.float64)[:, 1] - rcm[1]) * velocity.astype(np.float64)[:, 0])
                  )

        return angmom.sum()

    def max(self, quantity):
        return np.max(self.get_from_h5(quantity, save_to_mem=False))

    def min(self, quantity):
        return np.min(self.get_from_h5(quantity, save_to_mem=False))

    def num_particles(self, ptype=None):
        p = self.get_from_h5(h.NUMPARTTHISFILE)

        if ptype is None:
            return p.sum()
        else:
            return p[ptype]

    def coords(self, ptype=0):
        return self.get_from_h5(n.COORDINATES, ptype)

    def passive(self, ptype=0):
        return self.get_from_h5(n.PASSIVESCALARS, ptype)

    def vel(self, ptype=0):
        return self.get_from_h5(n.VELOCITIES, ptype)

    def mass(self, ptype=0):
        return self.get_from_h5(n.MASSES, ptype)

    def u(self, ptype=0):
        return self.get_from_h5(n.INTERNALENERGY, ptype)

    def xnuc(self, ptype=0):
        return self.get_from_h5(n.NUCLEARCOMPOSITION, ptype)

    def rho(self, ptype=0):
        return self.get_from_h5(n.DENSITY, ptype)

    def pressure(self, ptype=0):
        return self.get_from_h5(n.PRESSURE, ptype)

    def get_from_h5(self, field, ptype=0, from_file=False, save_to_mem=False):

        # Try to get from memory rather than disk first
        if not from_file:
            if field in self.data:
                # print(f"Getting {field} from self.data[field]")
                return self.data[field]

        p, isattr = path(field, ptype)

        with h5py.File(self.filename, 'r') as f:
            # print(f"Getting {field} from {self.filename}")

            if isattr:
                shape = f[p].attrs[field].shape

                if shape == ():
                    val = f[p].attrs[(field)]
                else:
                    val = np.array(f[p].attrs[field])
            else:
                val = np.array(f[p][field])

            if save_to_mem:
                self.data[field] = val
            return val

    def quick_pcolor(self, quantity_name, aa):

        boxsize = self.get_from_h5(n.BOXSIZE)
        t = self.get_from_h5(n.TIME)
        from run import ArepoRun
        ar = ArepoRun([self])

        if aa.inner_boxsize is not None:
            boxsize_x = aa.inner_boxsize
            boxsize_y = aa.inner_boxsize
            boxsize_z = aa.inner_boxsize
        else:
            boxsize_x = boxsize
            boxsize_y = boxsize
            boxsize_z = boxsize

        xlim = [boxsize / 2 - boxsize_x / 2, boxsize / 2 + boxsize_x / 2]
        ylim = [boxsize / 2 - boxsize_y / 2, boxsize / 2 + boxsize_y / 2]
        zlim = [boxsize / 2 - boxsize_z / 2, boxsize / 2 + boxsize_z / 2]


        cmap, norm, scmp = utilities.get_cmap(quantity_name, aa.cbar_lims, log_cmap=aa.log_cmap)

        units = part_fields[quantity_name]["Units"]

        po = PColorPlotOptions(
            quantity=quantity_name,
            orientation=["x", "y"],
            t=0,
            xlim=xlim,
            ylim=ylim,
            zlim=zlim,
            aa=aa,
            ar=ar,
            scmp=scmp,
            norm=norm,
            cmap=cmap,
            xlabel=f"x [cm]",
            ylabel=f"y [cm]",
            title=f"{quantity_name} [${units}$] (t={round(t, 2)})",
        )
        return PColorPlot(po)

    # def quick_pcolor_xnuc(self, quantity_name, export_filename=None, inner_boxsize=None, select_column=None):
    #
    #     res = 1024
    #     orien = ["x", "y"]
    #
    #     coords = self.get_from_h5(n.COORDINATES)
    #
    #     if select_column is not None:
    #         quant = self.get_from_h5(quantity_name)[:, select_column]
    #     else:
    #         quant = self.get_from_h5(quantity_name)
    #
    #     boxsize = self.get_from_h5(n.BOXSIZE)
    #
    #     if inner_boxsize is None:
    #         boxsize_x = boxsize
    #         boxsize_y = boxsize
    #         boxsizes = [boxsize_x, boxsize_y]
    #     else:
    #         boxsize_x = inner_boxsize
    #         boxsize_y = inner_boxsize
    #         boxsizes = [boxsize_x, boxsize_y]
    #
    #     resolutions = [res, res]
    #     xc, yc, _ = C.coordinates(orien)
    #     axes = [xc, yc]
    #     centers = [np.average([0, boxsize]), np.average([0, boxsize]), np.average([0, boxsize])]
    #
    #     with suppress_stdout_stderr():
    #         data = make_pcolor(coords.astype('float64'), quant.astype('float64'),
    #                            axes,
    #                            boxsizes,
    #                            resolutions,
    #                            centers,
    #                            include_neighbours_in_output=1,
    #                            numthreads=1)
    #
    #     x = np.arange(res + 1, dtype="float64") / res * boxsize_x - 0.5 * boxsize_x + centers[0]
    #     y = np.arange(res + 1, dtype="float64") / res * boxsize_y - 0.5 * boxsize_y + centers[1]
    #
    #     fig, ax = plt.subplots()
    #     pcolor = ax.pcolormesh(x, y, np.transpose(data["grid"]), vmin=0, vmax=1,
    #                            shading='flat', rasterized=True, cmap="hot")
    #     fig.colorbar(pcolor)
    #     ax.set(xlabel="x [cm]", ylabel="y [cm]", title=f"{quantity_name} (Carbon)")
    #
    #     if export_filename is not None:
    #         dir_name = os.path.dirname(self.filename)
    #         fn = os.path.join(dir_name, export_filename)
    #         fig.savefig(fn, dpi=300)
    #     else:
    #         plt.show()

    def quick_radial(self, quantity_name, aa):

        boxsize = self.get_from_h5(n.BOXSIZE)
        t = self.get_from_h5(n.TIME)
        from run import ArepoRun
        ar = ArepoRun([self])
        cmap, norm, scmp = utilities.get_cmap(quantity_name, aa.cbar_lims, log_cmap=aa.log_cmap)

        if aa.inner_boxsize is not None:
            boxsize_x = aa.inner_boxsize
            boxsize_y = aa.inner_boxsize
            boxsize_z = aa.inner_boxsize
        else:
            boxsize_x = boxsize
            boxsize_y = boxsize
            boxsize_z = boxsize

        xlim = [boxsize / 2 - boxsize_x / 2, boxsize / 2 + boxsize_x / 2]
        ylim = [boxsize / 2 - boxsize_y / 2, boxsize / 2 + boxsize_y / 2]
        zlim = [boxsize / 2 - boxsize_z / 2, boxsize / 2 + boxsize_z / 2]

        units = part_fields[quantity_name]["Units"]

        po = RadialPlotOptions(
            quantity=quantity_name,
            orientation=["x", "y"],
            t=0,
            xlim=xlim,
            ylim=ylim,
            zlim=zlim,
            aa=aa,
            ar=ar,
            scmp=scmp,
            norm=norm,
            cmap=cmap,
            xlabel=f"x [cm]",
            ylabel=f"y [cm]",
            title=f"{quantity_name} [${units}$] (t={round(t, 2)})",
            logscale=aa.logscale,
        )

        return RadialPlot(po)

    def quick_nuclear_compositions(self, export_filename=None):

        radials = []
        xnuc_names = ["Helium", "Carbon", "Oxygen", "Neon", "Magnesium"]
        colours = ["b", "k", "r", "g", "grey"]
        coords = self.get_from_h5(n.COORDINATES)
        quant = self.get_from_h5(n.NUCLEARCOMPOSITION)
        boxsize = self.get_from_h5(n.BOXSIZE)
        inner_boxsize = boxsize
        nspecies = np.shape(quant)[1]

        for x in range(nspecies):
            this_species = quant[:, x]

            y_av = 0.5 * boxsize
            z_av = 0.5 * boxsize
            boxsize_x = boxsize
            x_left = 0.5 * boxsize - 0.5 * inner_boxsize
            x_right = 0.5 * boxsize + 0.5 * inner_boxsize
            a = [x_left, y_av, z_av]
            b = [x_right, y_av, z_av]
            cyl_rad = 0.01 * boxsize
            nshells = 200
            radial = make_radial(coords.astype('float64'),
                                 this_species.astype('float64'),
                                 a, b,
                                 cyl_rad,
                                 nshells)

            radials.append(radial)

        # Make plot
        fig, ax = plt.subplots()
        lines = []

        for i, radial in enumerate(radials):
            line = ax.plot(radial[1, :], radial[0, :], colours[i], label=f"{xnuc_names[i]}")
            lines.append(line)

        the_sum = radials[0][0, :] + radials[1][0, :] + radials[2][0, :]
        ax.plot(radials[0][1, :], the_sum, 'yellow', label=f"sum")
        plt.ylim((0, 1.1))
        plt.xlabel("Radial distance [cm]")
        plt.ylabel(f"Nuclear Composition")
        plt.legend()
        plt.show()

        if export_filename is not None:
            dir_name = os.path.dirname(self.filename)
            fn = os.path.join(dir_name, export_filename)
            fig.savefig(fn, dpi=300)


class ArepoICs(ArepoH5File):

    def __init__(self, filename):
        super(ArepoICs, self).__init__(filename)

    def construct_header(self, particle_dict, config, param):

        num_part = np.zeros(6, dtype=np.int32)
        num_part[0] = np.shape(particle_dict[n.MASSES])[0]

        header = {
            h.NUMPARTTHISFILE: num_part,
            h.NUMPARTTOTAL: num_part,
            h.NUMPARTTOTALHIGHWORD: np.zeros(6, dtype=np.int32),
            h.MASSTABLE: np.zeros(6, dtype=np.float64),
            h.TIME: param.get("TimeBegin")["value"],
            h.BOXSIZE: param.get("BoxSize")["value"],
            h.REDSHIFT: np.array(0, dtype=np.float64),
            h.OMEGA0: np.array(param.get("Omega0")["value"], dtype=np.float64),
            h.OMEGALAMBDA: np.array(param.get("OmegaLambda")["value"], dtype=np.float64),
            h.HUBBLEPARAM: np.array(0, dtype=np.float64),
            h.NUMFILESPERSNAPSHOT: np.array(param.get("NumFilesPerSnapshot")["value"], dtype=np.int32),
            h.FLAGDOUBLEPRECISION: np.array(config.get("INPUT_IN_DOUBLEPRECISION")["value"], dtype=np.int32),
            h.FLAGFEEDBACK: np.array(0, dtype=np.int32),
            h.FLAGSFR: np.array(param.get("StarformationOn")["value"], dtype=np.int32),
            h.FLAGCOOLING: np.array(param.get("CoolingOn")["value"], dtype=np.int32),
            h.FLAGSTELLARAGE: np.array(0, dtype=np.int32),
            h.FLAGMETALS: np.array(0, dtype=np.int32),
            h.FLAGICINFO: np.array(0, dtype=np.int32),
            h.FLAGENTROPYICS: np.array(0, dtype=np.int32),
            }

        return header

    # config and param objects from Config and Param constructors
    def write_ics(self, particle_dict, config, param):

        print(f"Writing HDF5 file to: {self.filename}")

        # Construct the header
        header = self.construct_header(particle_dict, config, param)
        double_prec = bool(header['Flag_DoublePrecision'])

        with h5py.File(self.filename, "w") as f:

            # Write the header
            print(f"Writing header.")
            h = f.create_group(f"/{n.HEADER}")
            for key in header.keys():
                print(f"Writing key {key} with value {header[key]}")
                h.attrs[key] = header[key]

            i = 0
            print(f"Writing PartType{i}")
            name = f"/PartType{i}"
            nparticles = header[n.NUMPARTTOTAL][i]
            group = f.create_group(name)

            for quantity_name in particle_dict:
                dim = part_fields[quantity_name]["Dim"]
                quantity = particle_dict[quantity_name]

                # If dim is None, this means the dimensionality is variable. Check quantity.shape
                if dim is None:
                    dim = quantity.shape[1]

                if "ID" in quantity_name:
                    if double_prec:
                        this_dtype = np.uint64
                    else:
                        this_dtype = np.uint32
                else:
                    if double_prec:
                        this_dtype = np.float64
                    else:
                        this_dtype = np.float32

                print(f"Writing PartType{i}: {quantity_name}")
                group.create_dataset(quantity_name, shape=(nparticles, dim), dtype=this_dtype, data=quantity)

        print("Done.")

    @classmethod
    def from_snapshot(cls, snapshot_obj):
        raise NotImplementedError


class ArepoSnapshot(ArepoH5File):

    def __init__(self, filename):
        # self.data = dict()
        super(ArepoSnapshot, self).__init__(filename)

    def __str__(self):
        return f"ArepoSnapshot with n particles: {self.num_particles()}"
