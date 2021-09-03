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

        :return: The mean number of nucleons per isotope, the mean charge per isotope
        :rtype: (float, float)

        .. seealso::
            https://iopscience.iop.org/article/10.1086/313304/fulltext/40612.text.html
        .. notes::
            Let isotope i have Z protons and A nucleons (protons + neutrons).
            Let the aggregate total of isotope i have a number density of ni
            Define the dimensionless mass fraction of isotope i as:

            Xi = rhoi / rho = ni Ai / (rho NA)

            Then
            abar = Sum(Xi / Ai) ^(-1)
            zbar = abar * Sum(Zi Xi / Ai)

            Also:
            Ye = zbar / abar

        """
        abar_inv = 0.0
        zbar_sum = 0.0

        for i, spec in enumerate(self.species_dict):
            ai = self.species_dict[spec].mass_number
            zi = self.species_dict[spec].atomic_number
            xi = xnuc[i]

            abar_inv += xi / ai
            zbar_sum += zi * xi / ai

        abar = abar_inv ** (-1)
        zbar = abar * zbar_sum

        return abar, zbar

    def estimate_temp_from_e(self, e, xnuc, gamma=1.667):
        """Estimates temperature based on spec. internal energy and nuclear composition

        :param e: Specific internal energy
        :type e: float
        :param xnuc: Nuclear composition
        :type xnuc: np.ndarray
        :param gamma:
        :type gamma: float

        :return: Estimated temperature
        :rtype: float

        .. notes::
            From the C++ code
            _temp = 2.0 / 3.0 * p / (ni + ne) / GSL_CONST_CGS_BOLTZMANN;

            where
            - ni = 1.0 / cache.abar * rho * GSL_CONST_NUM_AVOGADRO;
            - ne = cache.zbar * ni;

            so:
            T = 2.0 / 3.0 * p / (ni + ne) / KB
            ni = rho * NA / abar
            ne = rho * NA / abar * zbar

            and we will use the relation
            p = e * rho * (gamma - 1)

            Thus
            T = 2.0 / 3.0 * e * rho * (gamma - 1) / (rho * NA / abar * (zbar + 1)) * 1 / KB
            T = 2.0 / 3.0 * e * (gamma - 1) / (NA * KB / abar * (zbar + 1))

            Rearranging for clarity:
            T = 2.0 / 3.0 * (abar * e * (gamma - 1)) / (KB * NA * (zbar + 1))
        """

        abar, zbar = self.azbar(xnuc)
        return 2.0 / 3.0 * (abar * e * (gamma - 1)) / (KB * NA * (zbar + 1))

    def estimate_e_from_temp(self, temp, xnuc, gamma=1.667):
        abar, zbar = self.azbar(xnuc)
        return 3.0 / 2.0 * (NA * KB * temp * (zbar + 1)) / (abar * (gamma - 1))
