from .names import ArepoHeader as h, n, path, ArepoGasFields, ArepoRoot
from .utilities import Coordinates as Coords, sci
from .sim_config import Param, Config
from .species import ArepoSpeciesList
from typing import Dict, Optional
import numpy as np
import h5py
import os


class ArepoH5File:
    """Base class for ArepoSnapshot and ArepoICs"""

    def __init__(self,
                 filename: str) -> None:
        """Constructor for hdf5 file object

        :param filename: Location of file on disk
        :type filename: str

        :return: ArepoH5File object
        :rtype: ArepoH5File
        """

        self.data = dict()
        self.maxima = dict()
        self.minima = dict()
        self.filename = os.path.abspath(filename)

    def r(self,
          center: np.ndarray) -> np.ndarray:
        """Returns an array of radii between center and each particle in Coordinates

        :param center: 3D coordinates of center
        :type center: np.ndarray

        :return: Array of radii
        :rtype: np.ndarray
        """
        coords = self.get_from_h5(n.COORDINATES)

        return np.sqrt(
            (coords[:, 0] - center[0]) ** 2 +
            (coords[:, 1] - center[1]) ** 2 +
            (coords[:, 2] - center[2]) ** 2)

    def center_of_mass(self) -> np.ndarray:
        """Computes the center of mass of all particles in Coordinates

        :return: Center of mass
        :rtype: np.ndarray
        """
        coords = self.get_from_h5(n.COORDINATES)
        mass = self.get_from_h5(n.MASSES)

        rcm = [(coords[:, i] / mass.sum() * mass).sum() for i in [0, 1, 2]]
        return np.array(rcm)

    def mean_a(self,
               species_list: ArepoSpeciesList) -> np.ndarray:
        """Computes mean molecular weight

        :param species_list:
        :type species_list: ArepoSpeciesList

        :return: Mean molecular weight array
        :rtype: np.ndarray
        """

        xnuc = self.get_from_h5(n.NUCLEARCOMPOSITION)

        mean_a = np.zeros(xnuc.shape[0])
        for i in range(species_list.num_species):
            a = list(species_list.species_dict.values())
            mean_a += xnuc[:, i] * a[i].mass_number
        return mean_a

    def coords_center_at(self,
                         center: np.ndarray) -> np.ndarray:
        """Returns an array which gives the offset from center to each particle in coordinates

        :param center: Center from which to calculate offsets
        :type center: np.ndarray

        :return: Array of offsets
        :rtype: np.ndarray
        """
        return self.get_from_h5(n.COORDINATES) - np.array(center)[None, :]

    def angular_momentum(self,
                         orientation: list[str]) -> np.ndarray:
        """Computes angular momentum given an orientation

        :param orientation: Orientation
        :type orientation: list[str]

        :return: Array of angular momentum values
        :rtype: np.ndarray
        """

        # Convert orientation into number representation
        xc, yc, _ = Coords.coordinates(orientation)

        # Get data
        velocity    = self.get_from_h5(n.VELOCITIES).astype(np.float64)
        mass        = self.get_from_h5(n.MASSES).astype(np.float64)
        rcm         = self.center_of_mass().astype(np.float64)
        coords      = self.coords_center_at(rcm).astype(np.float64)

        # Compute angular momentum
        angmom = mass * (+ (coords[:, xc] - rcm[xc]) * velocity[:, yc] - (coords[:, yc] - rcm[yc]) * velocity[:, xc])
        angmom = angmom.sum()

        return angmom

    def max(self,
            quantity: str) -> float:
        """Gets max value of quantity in the current file

        :param quantity: Quantity
        :type quantity: str

        :return: Max value of quantity
        :rtype: float
        """
        return np.max(self.get_from_h5(quantity, save_to_mem=False))

    def min(self,
            quantity: str) -> float:
        """Gets min value of quantity in the current file

        :param quantity: Quantity
        :type quantity: str

        :return: Min value of quantity
        :rtype: float
        """
        return np.min(self.get_from_h5(quantity, save_to_mem=False))

    def num_particles(self,
                      ptype: Optional[int] = None) -> int:
        """Gets the total number of particles in the file

        :param ptype: Particle type
        :type ptype: int

        :return: Number of particles of type ptype
        :rtype: int
        """

        p = self.get_from_h5(h.NUMPARTTHISFILE)

        if ptype is None:
            return p.sum()
        else:
            return p[ptype]

    def coords(self,
               ptype: int = 0) -> np.ndarray:
        """Convenience function to get coordinates.

        :param ptype: Particle type
        :type ptype: int

        :return: Coordinates
        :rtype: np.ndarray
        """
        return self.get_from_h5(n.COORDINATES, ptype)

    def passive(self,
                ptype: int = 0) -> np.ndarray:
        """Convenience function to get passive scalars.

        :param ptype: Particle type
        :type ptype: int

        :return: Passive scalars
        :rtype: np.ndarray
        """
        return self.get_from_h5(n.PASSIVESCALARS, ptype)

    def vel(self,
            ptype: int = 0) -> np.ndarray:
        """Convenience function to get velocities.

        :param ptype: Particle type
        :type ptype: int

        :return: Velocities
        :rtype: np.ndarray
        """
        return self.get_from_h5(n.VELOCITIES, ptype)

    def mass(self,
             ptype: int = 0) -> np.ndarray:
        """Convenience function to get masses.

        :param ptype: Particle type
        :type ptype: int

        :return: Masses
        :rtype: np.ndarray
        """
        return self.get_from_h5(n.MASSES, ptype)

    def u(self,
          ptype: int = 0) -> np.ndarray:
        """Convenience function to get internal energy.

        :param ptype: Particle type
        :type ptype: int

        :return: Internal energy
        :rtype: np.ndarray
        """
        return self.get_from_h5(n.INTERNALENERGY, ptype)

    def xnuc(self,
             ptype: int = 0) -> np.ndarray:
        """Convenience function to get nuclear composition.

        :param ptype: Particle type
        :type ptype: int

        :return: Nuclear composition
        :rtype: np.ndarray
        """
        return self.get_from_h5(n.NUCLEARCOMPOSITION, ptype)

    def rho(self,
            ptype: int = 0) -> np.ndarray:
        """Convenience function to get density.

        :param ptype: Particle type
        :type ptype: int

        :return: Density
        :rtype: np.ndarray
        """
        return self.get_from_h5(n.DENSITY, ptype)

    def pressure(self,
                 ptype: int = 0) -> np.ndarray:
        """Convenience function to get pressure.

        :param ptype: Particle type
        :type ptype: int

        :return: Pressure
        :rtype: np.ndarray
        """
        return self.get_from_h5(n.PRESSURE, ptype)

    def get_from_h5(self,
                    field: str,
                    ptype: int = 0,
                    from_file: bool = False,
                    save_to_mem: bool = False):
        """General function to get some data from the h5 file or the internal self.data cache

        :param field: Field to obtain data from
        :type field: str
        :param ptype: Particle type
        :type ptype: int
        :param from_file: Source of data, either from file or cache
        :type from_file: bool
        :param save_to_mem: If data should be saved to internal cache
        :type save_to_mem: bool

        :return: Data from the h5 file
        :rtype: np.ndarray
        """

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

    def particle_exists(self,
                        pid: int) -> bool:
        pids = self.get_from_h5(n.PARTICLEIDS)
        return np.isin(pids, np.array([pid])).any()

    def get_ids_of_largest(self,
                           quant_name: str,
                           k: int = 10,
                           pids: list[int] = None) -> list[int]:

        if pids is None:
            pids = self.get_from_h5(n.PARTICLEIDS)

        quant = self.get_values_from_pids(quant_name, pids)
        return self.get_pids_of_k_largest_from_array(pids, quant, k)

    def get_ids_of_smallest(self,
                            quant_name: str,
                            k: int = 10,
                            pids: list[int] = None) -> list[int]:

        if pids is None:
            pids = self.get_from_h5(n.PARTICLEIDS)

        quant = self.get_values_from_pids(quant_name, pids)
        return self.get_pids_of_k_smallest_from_array(pids, quant, k)

    def get_values_from_pids(self,
                             quant_name: str,
                             pids_requested: np.ndarray) -> np.ndarray:

        quant = self.get_from_h5(quant_name)
        pid = self.get_from_h5(n.PARTICLEIDS)
        mask = np.isin(pid, pids_requested)

        # Error if non-existant id requested
        if np.count_nonzero(mask) != len(pids_requested):
            found_particles = pid[mask]
            not_found_particles = np.setdiff1d(pids_requested, found_particles)
            raise ValueError(f"Found particles {found_particles} but not {not_found_particles}")

        orig_indices = pid.argsort()
        ndx = orig_indices[np.searchsorted(pid[orig_indices], pids_requested)]

        return quant[ndx]

    def get_field_names(self) -> list[str]:
        field_names = []
        for field in ArepoGasFields.values():
            try:
                self.get_from_h5(field)
                field_names.append(field)
            except (KeyError, ValueError):
                pass

        return field_names

    def print_particle_info(self,
                            list_of_pids: np.ndarray,
                            limit: int = 10,
                            fields: list[str] = None) -> None:

        # Limit the output
        if len(list_of_pids) > limit:
            list_of_pids = list_of_pids[:limit]

        if fields is None:
            fields = self.get_field_names()

        info = np.full((len(list_of_pids), ), "", dtype=f"S{50 * len(fields)}")

        for quant_name in fields:
            vals = self.get_values_from_pids(quant_name, list_of_pids)
            for i, p in enumerate(list_of_pids):
                if vals.ndim == 1:
                    formatted_val = sci(vals[i])
                else:
                    formatted_val = vals[i, :]

                info[i] += str.encode(f"{quant_name}: {formatted_val} \n")

        np.set_printoptions(precision=4)
        for i, info_str in enumerate(info):
            print(f"{list_of_pids[i]=}")
            print(info_str.decode())

    @staticmethod
    def get_pids_of_k_largest_from_array(pids: list[int],
                                         quant: np.ndarray,
                                         k: int):
        return pids[np.argpartition(quant, -k)[-k:]]

    @staticmethod
    def get_pids_of_k_smallest_from_array(pids: list[int],
                                          quant: np.ndarray,
                                          k: int):
        return pids[np.argpartition(quant, k)[:k]]


class ArepoICs(ArepoH5File):
    """Class to hold information for Arepo initial conditions."""

    def __init__(self,
                 filename: str) -> None:
        """Constructor for ArepoICs

        :param filename: Filename of IC file
        :type filename: str

        :return: ArepoICs object
        :rtype: ArepoICs
        """
        super(ArepoICs, self).__init__(filename)

    def write_ics(self,
                  particle_dict: dict,
                  config: Config,
                  param: Param,
                  add_pids: bool = True,
                  long_ids: bool = False) -> None:
        """Write the initial conditions information to file.

        :param particle_dict: Particle dictionary which specifies coordinates, velocities, etc.
        :type particle_dict: dict
        :param config: Config object
        :type config: Config
        :param param: Param object
        :type param: Param
        :param add_pids: Switch to add ParticleIDs to file
        :type add_pids: bool
        :param long_ids: Switch to write ParticleIDs as np.uint64
        :type long_ids: bool
        """

        print(f"Writing HDF5 file to: {self.filename}")

        # Construct the header
        header_dict = self.construct_header(particle_dict, config, param)
        double_prec = bool(header_dict[h.FLAGDOUBLEPRECISION])

        with h5py.File(self.filename, "w") as f:

            # Write the header
            print("Writing header.")
            header = f.create_group(f"/{ArepoRoot.HEADER}")
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
    def construct_header(particle_dict: dict,
                         config: Config,
                         param: Param) -> Dict[str, np.ndarray]:
        """Construct a header for writing IC to file.

        :param particle_dict: Particle dictionary
        :type particle_dict: dict
        :param config: Config object
        :type config: Config
        :param param: Param object
        :type param: Param

        :return: Header
        :rtype: dict
        """

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
    """ArepoSnapshot."""

    def __init__(self,
                 filename: str) -> None:
        """Constructor for Snapshots

        :param filename: Filename of snapshot
        :type filename: str

        :return: Areposnapshot object
        :rtype: ArepoSnapshot
        """
        super(ArepoSnapshot, self).__init__(filename)

    def __str__(self) -> str:
        return f"ArepoSnapshot with n particles: {self.num_particles()}"
