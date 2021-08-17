from names import ArepoHeader as h, n, path
from utilities import part_fields
import numpy as np
import h5py
import os


class ArepoH5File(object):

    def __init__(self, filename):
        self.data = dict()
        self.maxima = dict()
        self.minima = dict()
        self.filename = os.path.abspath(filename)

    def r(self, center):
        coords = self.get_from_h5(n.COORDINATES)

        return np.sqrt(
            (coords[:, 0] - center[0]) ** 2 +
            (coords[:, 1] - center[1]) ** 2 +
            (coords[:, 2] - center[2]) ** 2)

    def center_of_mass(self):
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
                return self.data[field]

        p, isattr = path(field, ptype)

        with h5py.File(self.filename, "r") as f:

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


class ArepoICs(ArepoH5File):

    def __init__(self, filename):
        super(ArepoICs, self).__init__(filename)

    def write_ics(self, particle_dict, config, param, add_pids=True, long_ids=False):

        print(f"Writing HDF5 file to: {self.filename}")

        # Construct the header
        header_dict = self.construct_header(particle_dict, config, param)
        double_prec = bool(header_dict[h.FLAGDOUBLEPRECISION])

        with h5py.File(self.filename, "w") as f:

            # Write the header
            print(f"Writing header.")
            header = f.create_group(f"/{n.HEADER}")
            for key in header_dict.keys():
                print(f"Writing key {key} with value {header_dict[key]}")
                header.attrs[key] = header_dict[key]

            i = 0
            print(f"Writing PartType{i}")
            name = f"/PartType{i}"
            nparticles = header_dict[n.NUMPARTTOTAL][i]
            group = f.create_group(name)

            if add_pids:
                particle_dict[n.PARTICLEIDS] = np.arange(1, nparticles + 1)

            for quantity_name in particle_dict:
                # dim = part_fields[quantity_name]["Dim"]
                quantity = particle_dict[quantity_name]

                # If dim is None, this means the dimensionality is variable. Check quantity.shape
                if len(quantity.shape) > 1:
                    dim = quantity.shape[1]
                    shape = (nparticles, dim)
                else:
                    shape = (nparticles, )

                if "ID" in quantity_name:
                    if long_ids:
                        this_dtype = np.uint64
                    else:
                        this_dtype = np.uint32
                else:
                    if double_prec:
                        this_dtype = np.float64
                    else:
                        this_dtype = np.float32

                print(f"Writing PartType{i}: {quantity_name}")
                group.create_dataset(quantity_name, shape=shape, dtype=this_dtype, data=quantity)

        print("Done.")

    @staticmethod
    def construct_header(particle_dict, config, param):

        num_part = np.zeros(6, dtype=np.int32)
        num_part[0] = np.shape(particle_dict[n.MASSES])[0]

        header = {
            h.NUMPARTTHISFILE: num_part,
            h.NUMPARTTOTAL: num_part,
            h.NUMPARTTOTALHIGHWORD: np.zeros(6, dtype=np.int32),
            h.MASSTABLE: np.zeros(6, dtype=np.float64),
            h.TIME: param.get("TimeBegin"),
            h.BOXSIZE: param.get(h.BOXSIZE),
            h.REDSHIFT: np.array(0, dtype=np.float64),
            h.OMEGA0: np.array(param.get(h.OMEGA0), dtype=np.float64),
            h.OMEGALAMBDA: np.array(param.get(h.OMEGALAMBDA), dtype=np.float64),
            h.HUBBLEPARAM: np.array(0, dtype=np.float64),
            h.NUMFILESPERSNAPSHOT: np.array(param.get(h.NUMFILESPERSNAPSHOT), dtype=np.int32),
            h.FLAGDOUBLEPRECISION: np.array(config.get("INPUT_IN_DOUBLEPRECISION"), dtype=np.int32),
            h.FLAGFEEDBACK: np.array(0, dtype=np.int32),
            h.FLAGSFR: np.array(param.get("StarformationOn"), dtype=np.int32),
            h.FLAGCOOLING: np.array(param.get("CoolingOn"), dtype=np.int32),
            h.FLAGSTELLARAGE: np.array(0, dtype=np.int32),
            h.FLAGMETALS: np.array(0, dtype=np.int32),
            h.FLAGICINFO: np.array(0, dtype=np.int32),
            h.FLAGENTROPYICS: np.array(0, dtype=np.int32),
            }

        return header

    @classmethod
    def from_snapshot(cls, snapshot_obj):
        raise NotImplementedError


class ArepoSnapshot(ArepoH5File):

    def __init__(self, filename):
        super(ArepoSnapshot, self).__init__(filename)

    def __str__(self):
        return f"ArepoSnapshot with n particles: {self.num_particles()}"
