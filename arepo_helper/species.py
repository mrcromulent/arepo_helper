from collections import OrderedDict
from mendeleev import element
from const import NA, KB
import numpy as np
import os


class ArepoSpeciesList(object):
    """Class to work with species.txt files for the nuclear network"""

    species_dict    = OrderedDict()
    num_species     = 0

    def __init__(self, species_file):
        """Constructor to load data from the species.txt file.

        :param species_file: Species file name
        :type species_file: str

        :raises FileNotFoundError: If species file not found
        :rtype: None
        """

        if not os.path.exists(species_file):
            raise FileNotFoundError(f"Species file {species_file} does not exist.")

        self.load_species(species_file)

    def index_of(self, elem):
        """Index of a particular species

        :param elem: Element name (e.g. He4)
        :type elem: str

        :return: Index of the element in the file
        :rtype: int
        """
        return list(self.species_dict).index(elem)

    def load_species(self, filename):
        """Loads the data located in the file at filename

        :param filename: species.txt filename
        :type filename: str
        """

        data = OrderedDict()
        with open(filename) as f:

            # Skip the header row
            next(f)
            for line in f.readlines():
                if not line == "\n":
                    name, na, nz = line.split()
                    name    = name.capitalize()
                    na      = int(na)
                    name_no_na = "".join(i for i in name if not i.isdigit())
                    elem = element(name_no_na)

                    selected_isotope = None
                    for iso in elem.isotopes:
                        if iso.mass_number == na:
                            selected_isotope = iso

                    data[name] = selected_isotope

        self.species_dict   = data
        self.num_species    = len(data.keys())

    def azbar(self, xnuc):
        """Computes the mean number of baryons and mean charge per baryon

        :param xnuc: Nuclear composition, indicating the percentage of each species
        :type xnuc: np.ndarray
        :return: Mean baryons and mean charge per baryon
        :rtype: (float, float)
        """
        abar = 0.0
        zbar = 0.0
        xsum = 0.0

        for i, spec in enumerate(self.species_dict):
            curr_a = self.species_dict[spec].mass_number
            curr_z = self.species_dict[spec].atomic_number

            ymass = xnuc[i] * curr_a
            abar += ymass
            zbar += curr_z * ymass
            xsum += xnuc[i]

        abar = xsum / abar
        zbar = zbar / xsum * abar

        return abar, zbar

    def estimate_temp(self, rho, e, xnuc):
        """Estimates temperature based on density, spec. internal energy and nuclear composition

        :param rho: Density
        :type rho: float
        :param e: Specific internal energy
        :type e: float
        :param xnuc: Nuclear composition
        :type xnuc: np.ndarray

        :return: Estimated temperature
        :rtype: float
        """

        abar, zbar = self.azbar(xnuc)
        ni = 1 / abar * rho * NA
        ne = zbar * ni

        temp = 2.0 / (3.0 * (ni + ne)) * 1 / KB * e * rho

        return temp
