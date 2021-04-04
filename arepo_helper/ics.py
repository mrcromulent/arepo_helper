from utilities import header_default, part_fields
from createICs import create_particles_fill_grid
from h5_file import ArepoH5File
from names import n
import numpy as np
import h5py


class ArepoICs(ArepoH5File):

    def __init__(self, filename):
        super(ArepoICs, self).__init__(filename)

    def add_grids(self, grids, particle_dict, header, resolution, xnuc=None):
        for boxsize in grids:
            particle_dict, header = self.add_grid(particle_dict, header, boxsize, resolution, xnuc)

        return particle_dict, header

    def add_grid(self, pd, header, boxsize_to_add, resolution, xnuc):
        """
        :param pd:
        :param header:
        :param boxsize_to_add:
        :param resolution:
        :param xnuc:
        :return:

        Essentially the function embeds an existing particle dictionary inside a larger one, which is lower resolution.
        The larger region is populated with particles too but at a much lower density.
        """

        if n.BOXSIZE not in header:
            raise ValueError(f"{n.BOXSIZE} must be included in the header")

        p0 = pd[f"{n.PARTTYPE}0"]
        p0[n.COORDINATES] += 0.5 * boxsize_to_add - 0.5 * header[n.BOXSIZE]
        f64coords = p0[n.COORDINATES].astype('float64')

        # Update the number of particles
        p = create_particles_fill_grid(f64coords, boxsize_to_add, resolution)
        npart_old = header[n.NUMPARTTOTAL].sum()
        npart_new = npart_old + np.shape(p)[0]

        for q in [n.COORDINATES, n.VELOCITIES, n.MAGNETICFIELD, n.MASSES, n.INTERNALENERGY, n.NUCLEARCOMPOSITION,
                  n.PASSIVESCALARS, n.PARTICLEIDS]:
            if q in p0:
                sh = np.shape(p0[q])
                if len(sh) > 1:
                    p0[q] = np.resize(p0[q], (npart_new, sh[1]))
                    # p0[q].resize((npart_new, sh[1]))
                else:
                    p0[q] = np.resize(p0[q], npart_new)
                    # p0[q].resize(npart_new)

        p0[n.COORDINATES][npart_old:, :] = p

        if n.NUCLEARCOMPOSITION in p0:
            if xnuc is not None:
                p0[n.NUCLEARCOMPOSITION][npart_old:, :] = xnuc[None, :]

        header[n.NUMPARTTHISFILE] = np.array([npart_new, 0, 0, 0, 0, 0])
        header[n.NUMPARTTOTAL] = np.array([npart_new, 0, 0, 0, 0, 0])
        header[n.NUMPARTTOTALHIGHWORD] = np.array([npart_new, 0, 0, 0, 0, 0])

        return pd, header

    def write_ics(self, particle_dict, *args, **kwargs):
        """
        :param particle_dict: An ORDERED dictionary containing the important fields in the following format:
        {
        "PartType0": {"Coordinates": np.array(...), "Masses": np.array(...), ...},
        "PartType1": {"Coordinates": np.array(...), "Masses": np.array(...), ...},
        ...
        }

        args contains the header fields the user wishes to override from the default
        kwargs are used for additional information such as whether or not to use long uint ids
        """

        print(f"Writing HDF5 file to: {self.filename}")

        # Extra config info
        config = {
            "long_ids": False,
            "start_time": 0,
            "double": False}

        for key in kwargs:
            config[key] = kwargs[key]

        dtype = np.float64 if config["double"] else np.float32
        dtypeids = np.uint64 if config["long_ids"] else np.uint32

        # The header contains all the relevant information that goes into the HDF5 file
        header = header_default
        for dictionary in args:
            for key in dictionary:
                header[key] = dictionary[key]

        with h5py.File(self.filename, "w") as f:

            # Write the header
            print(f"Writing header.")
            h = f.create_group(f"/{n.HEADER}")
            for key in header.keys():
                print(f"Writing key {key} with value {header[key]}")
                h.attrs[key] = header[key]

            for i, part_type in enumerate(particle_dict.keys()):
                print(f"Writing PartType{i}")
                name = f"/PartType{i}"
                nparticles = header[n.NUMPARTTOTAL][i]
                group = f.create_group(name)

                for quantity_name in particle_dict[part_type]:
                    dim = part_fields[quantity_name]["Dim"]
                    quantity = particle_dict[part_type][quantity_name]

                    # If dim is None, this means the dimensionality is variable. Check quantity.shape
                    if dim is None:
                        dim = quantity.shape[1]

                    if "ID" in quantity_name:
                        this_dtype = dtypeids
                    else:
                        this_dtype = dtype

                    print(f"Writing PartType{i}: {quantity_name}")
                    group.create_dataset(quantity_name, shape=(nparticles, dim), dtype=this_dtype, data=quantity)

        print("Done.")

    @staticmethod
    def gadget_write_ics_format3(filename,
                                 data,
                                 time=0.0,
                                 double=False,
                                 longids=False,
                                 num_files=1,
                                 boxsize=0.0,
                                 masses=None,
                                 skipxnuc=False):
        import h5py

        filename += '.hdf5'
        print("Writing gadget file: ", filename)
        f = h5py.File(filename, 'w')

        npart = np.zeros(6, dtype=np.int32)
        nparthighword = np.zeros(6, dtype=np.int32)
        offset = np.zeros(6, dtype=np.int32)
        npartmass = np.zeros(6, dtype=np.int32)
        massoffset = np.zeros(6, dtype=np.int32)

        if 'type' in data:
            for ptype in range(6):
                npart[ptype] = np.size(np.where(data['type'] == ptype))
        else:
            npart[0] = data['count']

        offset[1:] = np.cumsum(npart[:-1])

        if not masses is None:
            massarr = masses
        else:
            massarr = np.zeros(6, dtype=np.float64)

        npartmass[:] = npart[:]
        j, = np.where(massarr > 0.0)
        npartmass[j] = 0
        massoffset[1:] = np.cumsum(npartmass[:-1])

        header = f.create_group("/Header")

        header.attrs['NumPart_ThisFile'] = npart
        header.attrs['NumPart_Total'] = npart
        header.attrs['NumPart_Total_HighWord'] = nparthighword
        header.attrs['MassTable'] = massarr
        header.attrs['Time'] = time
        header.attrs['NumFilesPerSnapshot'] = np.array(num_files, dtype=np.int32)
        header.attrs['BoxSize'] = boxsize
        header.attrs['Flag_DoublePrecision'] = np.array(double, dtype=np.int32)
        header.attrs['Flag_IC_info'] = np.array(0, dtype=np.int32)
        header.attrs['Flag_Entropy_ICs'] = np.array(0, dtype=np.int32)
        header.attrs['Redshift'] = np.array(0, dtype=np.float64)
        header.attrs['Omega0'] = np.array(0, dtype=np.float64)
        header.attrs['OmegaLambda'] = np.array(0, dtype=np.float64)
        header.attrs['HubbleParam'] = np.array(0, dtype=np.float64)
        header.attrs['Flag_Sfr'] = np.array(0, dtype=np.int32)
        header.attrs['Flag_Cooling'] = np.array(0, dtype=np.int32)
        header.attrs['Flag_StellarAge'] = np.array(0, dtype=np.int32)
        header.attrs['Flag_Metals'] = np.array(0, dtype=np.int32)
        header.attrs['Flag_Feedback'] = np.array(0, dtype=np.int32)

        if double:
            dtype = np.float64
        else:
            dtype = np.float32

        if longids:
            dtypeids = np.uint64
        else:
            dtypeids = np.uint32

        fields = []
        if 'bfld' in data:
            fields += ['bfld']
        if not skipxnuc and ('xnuc' in data):
            fields += ['xnuc']
        if 'temp' in data:
            fields += ['temp']
        if 'rho' in data:
            fields += ["rho"]
        if 'pass' in data:
            fields += ["pass"]
        if 'erad' in data:
            fields += ["erad"]

        fields_to_names = {
            'bfld': 'MagneticField',
            'xnuc': 'NuclearComposition',
            'temp': 'Temperature',
            'rho': 'Density',
            'pass': 'PassiveScalars',
            'erad': 'Erad'
        }

        for ptype in range(6):
            if npart[ptype] > 0:
                group = f.create_group('/PartType%d' % ptype)

                pos = group.create_dataset("Coordinates", (npart[ptype], 3), dtype)
                pos[:, :] = data['pos'][offset[ptype]:offset[ptype] + npart[ptype], :]

                vel = group.create_dataset("Velocities", (npart[ptype], 3), dtype)
                vel[:, :] = data['vel'][offset[ptype]:offset[ptype] + npart[ptype], :]

                id = group.create_dataset("ParticleIDs", (npart[ptype],), dtypeids)
                if 'id' in data:
                    id[:] = data['id'][offset[ptype]:offset[ptype] + npart[ptype]]
                else:
                    id[:] = np.arange(offset[ptype] + 1, offset[ptype] + npart[ptype] + 1, dtype=dtypeids)

                if massarr[ptype] == 0:
                    mass = group.create_dataset("Masses", (npart[ptype],), dtype)
                    mass[:] = data['mass'][massoffset[ptype]:massoffset[ptype] + npart[ptype]]

                if ptype == 0:
                    u = group.create_dataset("InternalEnergy", (npart[ptype],), dtype)
                    u[:] = data['u'][offset[ptype]:offset[ptype] + npart[ptype]]
                    # loop over the other fields
                    for field in fields:
                        hdf5key = fields_to_names[field]
                        if field in ['bfld', 'xnuc', 'pass']:
                            val = group.create_dataset(hdf5key, (npart[ptype], data[field].shape[1]), dtype)
                            val[:, :] = data[field][offset[ptype]:offset[ptype] + npart[ptype], :]
                        else:
                            val = group.create_dataset(hdf5key, (npart[ptype],), dtype)
                            val[:] = data[field][offset[ptype]:offset[ptype] + npart[ptype]]
        f.close()

        print("Done.")
        return

    @classmethod
    def from_snapshot(cls, snapshot_obj):
        pass


if __name__ == "__main__":
    test = ArepoH5File("./snapshots/snapshot_024.hdf5")

    pdict = {f"{n.PARTTYPE}0": {n.COORDINATES: np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, 5]])}}
    h = dict()
    h[n.NUMPARTTHISFILE] = np.array([2, 0, 0, 0, 0, 0])
    h[n.NUMPARTTOTAL] = np.array([2, 0, 0, 0, 0, 0])
    h[n.NUMPARTTOTALHIGHWORD] = np.array([2, 0, 0, 0, 0, 0])

    a = ArepoICs("./snapshots/test.hdf5")
    a.write_ics(pdict, h)
