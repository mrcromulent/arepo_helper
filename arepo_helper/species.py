from collections import OrderedDict
from mendeleev import element
import os


class ArepoSpeciesList(object):

    species_dict    = OrderedDict()
    num_species     = 0

    def __init__(self, species_file):

        if not os.path.exists(species_file):
            raise FileNotFoundError(f"Species file {species_file} does not exist.")

        self.load_species(species_file)

    def index_of(self, elem):
        return list(self.species_dict).index(elem)

    def load_species(self, filename):

        data = OrderedDict()
        with open(filename) as f:

            # Skip the header row
            next(f)

            for line in f.readlines():
                if not line == "\n":
                    name, na, nz = line.split()
                    name    = name.capitalize()
                    na      = int(na)
                    name_no_na = ''.join(i for i in name if not i.isdigit())
                    elem = element(name_no_na)

                    selected_isotope = None
                    for iso in elem.isotopes:
                        if iso.mass_number == na:
                            selected_isotope = iso

                    data[name] = selected_isotope

        self.species_dict   = data
        self.num_species    = len(data.keys())
