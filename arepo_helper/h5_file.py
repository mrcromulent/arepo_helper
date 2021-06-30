"""
Tried to use hdfdict but it flat out didn't work. Did not dump files correctly nor did it preserve attributes of groups
"""

from names import ArepoHeader, n, path
import numpy as np
import h5py
import os


class ArepoH5File(object):

    def __init__(self, filename):
        self.data       = dict()
        self.maxima     = dict()
        self.minima     = dict()
        self.filename   = os.path.abspath(filename)

    def r(self, center):
        # TODO: Adapt for lower dimensionality

        coords = self.get_from_h5(n.COORDINATES)

        return np.sqrt(
            (coords[:, 0] - center[0]) ** 2 +
            (coords[:, 1] - center[1]) ** 2 +
            (coords[:, 2] - center[2]) ** 2)

    def center_of_mass(self):
        # TODO: Adapt for lower dimensionality

        coords  = self.get_from_h5(n.COORDINATES)
        mass    = self.get_from_h5(n.MASSES)

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

        velocity    = self.get_from_h5(n.VELOCITIES)
        mass        = self.get_from_h5(n.MASSES)
        rcm         = self.center_of_mass()
        coords      = self.coords_center_at(rcm)

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
        p = self.get_from_h5(ArepoHeader.NUMPARTTHISFILE)

        if ptype is None:
            return p.sum()
        else:
            return p[ptype]

    def coords(self, ptype=0):
        return self.get_from_h5(n.COORDINATES, ptype)

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
