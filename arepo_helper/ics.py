from utilities import header_default, part_fields
from h5_file import ArepoH5File
from names import n
import numpy as np
import h5py


class ArepoICs(ArepoH5File):

    def __init__(self, filename):
        super(ArepoICs, self).__init__(filename)

    def add_grid(self, particle_dict, header, boxsize, resolution, xnuc):
        raise NotImplementedError

    def add_grids(self, grids, particle_dict, header, resolution, xnuc=None):
        for boxsize in grids:
            particle_dict, header = self.add_grid(particle_dict, header, boxsize, resolution, xnuc)

        return particle_dict, header

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

    @classmethod
    def from_snapshot(cls, snapshot_obj):
        raise NotImplementedError
