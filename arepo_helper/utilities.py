from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize, SymLogNorm
from names import ArepoHeader as h, n
import numpy as np
import os


class Coordinates:
    x = 0
    y = 1
    z = 2

    @classmethod
    def coordinates(cls, orientation):
        xstr, ystr = orientation
        xval = getattr(cls, xstr)
        yval = getattr(cls, ystr)
        zset = {cls.x, cls.y, cls.z} - {xval, yval}

        return xval, yval, list(zset)[0]


common_snapbases = ["snapshot", "snap", "sh"]
valid_orientations = [["x", "y"], ["x", "z"], ["y", "z"]]
default_cmap = "viridis"
plot_quantities = [n.DENSITY, n.TEMPERATURE, n.PRESSURE, n.MASSES, n.MAGNETICFIELD]
# plot_quantities = [n.DENSITY]

header_default = {h.NUMPARTTHISFILE: np.zeros(6, dtype=np.int32),
                  h.NUMPARTTOTAL: np.zeros(6, dtype=np.int32),
                  h.NUMPARTTOTALHIGHWORD: np.zeros(6, dtype=np.int32),
                  h.MASSTABLE: np.zeros(6, dtype=np.float64),
                  h.TIME: 0.,
                  h.NUMFILESPERSNAPSHOT: 1,
                  h.BOXSIZE: 0.,
                  h.FLAGDOUBLEPRECISION: np.array(False, dtype=np.int32),
                  h.FLAGICINFO: np.array(0, dtype=np.int32),
                  h.FLAGENTROPYICS: np.array(0, dtype=np.int32),
                  h.REDSHIFT: np.array(0, dtype=np.float64),
                  h.OMEGA0: np.array(0, dtype=np.float64),
                  h.OMEGALAMBDA: np.array(0, dtype=np.float64),
                  h.HUBBLEPARAM: np.array(0, dtype=np.float64),
                  h.FLAGSFR: np.array(0, dtype=np.int32),
                  h.FLAGSTELLARAGE: np.array(0, dtype=np.int32),
                  h.FLAGMETALS: np.array(0, dtype=np.int32),
                  h.FLAGFEEDBACK: np.array(0, dtype=np.int32),
                  h.FLAGCOOLING: np.array(0, dtype=np.int32)}

part_fields = {n.CENTEROFMASS: {"Dim": None, "Units": "cm", "cmap": default_cmap},
               n.COORDINATES: {"Dim": None, "Units": "cm", "cmap": default_cmap},
               n.DENSITY: {"Dim": 1, "Units": "g / cm^3", "cmap": default_cmap},
               n.DENSITYGRADIENT: {"Dim": None, "Units": "g / cm^3 / cm", "cmap": default_cmap},
               n.INTERNALENERGY: {"Dim": 1, "Units": "erg", "cmap": default_cmap},
               n.MAGNETICFIELD: {"Dim": None, "Units": "10^4 G", "cmap": default_cmap},
               n.MAGNETICFIELDDIVERGENCE: {"Dim": 1, "Units": "10^4 G / cm^3", "cmap": default_cmap},
               n.MASSES: {"Dim": 1, "Units": "g", "cmap": default_cmap},
               n.NUCLEARCOMPOSITION: {"Dim": None, "Units": None, "cmap": default_cmap},
               n.NUCLEARENERGYGENERATIONRATE: {"Dim": None, "Units": "erg", "cmap": default_cmap},
               n.PARTICLEIDS: {"Dim": 1, "Units": None, "cmap": default_cmap},
               n.PASSIVESCALARS: {"Dim": 4, "Units": None, "cmap": default_cmap},
               n.POTENTIAL: {"Dim": 1, "Units": "erg/g", "cmap": default_cmap},
               n.PRESSURE: {"Dim": 1, "Units": "dynes", "cmap": default_cmap},
               n.PRESSUREGRADIENT: {"Dim": None, "Units": "dynes/cm", "cmap": default_cmap},
               n.TEMPERATURE: {"Dim": 1, "Units": "K", "cmap": default_cmap},
               n.VELOCITIES: {"Dim": None, "Units": "cm/s", "cmap": default_cmap},
               n.VELOCITYDIVERGENCE: {"Dim": 1, "Units": "1/s/cm^2", "cmap": default_cmap}}


def get_cmap(quantity, bounds, log_cmap=False):
    """
    :param quantity: One of the particle fields from the enum
    :param bounds: A list containing the minimum value in the 0th position and the maximum in the first
    :return:
    """

    assert quantity in plot_quantities, f"Cmap requested for unknown quantity {quantity}"

    cmap = part_fields[quantity]["cmap"]

    if log_cmap:
        norm = SymLogNorm(10**5, vmin=bounds[0], vmax=bounds[1], base=10)
    else:
        norm = Normalize(bounds[0], bounds[1])
    scmp = ScalarMappable(norm=norm, cmap=cmap)

    return cmap, norm, scmp


def dummy():
    pass


class suppress_stdout_stderr(object):
    """
    A context manager for doing a "deep suppression" of stdout and stderr in
    Python, i.e. will suppress all print, even if the print originates in a
    compiled C/Fortran sub-function.
       This will not suppress raised exceptions, since exceptions are printed
    to stderr just before a script exits, and after the context manager has
    exited (at least, I think that is why it lets exceptions through).
    """

    def __init__(self):
        # Open a pair of null files
        self.null_fds = [os.open(os.devnull, os.O_RDWR) for x in range(2)]
        # Save the actual stdout (1) and stderr (2) file descriptors.
        self.save_fds = [os.dup(1), os.dup(2)]

    def __enter__(self):
        # Assign the null pointers to stdout and stderr.
        os.dup2(self.null_fds[0], 1)
        os.dup2(self.null_fds[1], 2)

    def __exit__(self, *_):
        # Re-assign the real stdout/stderr back to (1) and (2)
        os.dup2(self.save_fds[0], 1)
        os.dup2(self.save_fds[1], 2)
        # Close all file descriptors
        for fd in self.null_fds + self.save_fds:
            os.close(fd)


ruediger_conversion_table = {n.COORDINATES: "pos",
                             n.VELOCITIES: "vel",
                             n.NUCLEARCOMPOSITION: "xnuc",
                             n.INTERNALENERGY: "u",
                             n.PARTICLEIDS: "id",
                             n.MASSES: "mass",
                             n.VOLUME: "vol",
                             n.PRESSURE: "pres",
                             n.PASSIVESCALARS: "pass",
                             n.DENSITY: "rho",
                             n.TEMPERATURE: "temp",
                             n.POTENTIAL: "pot"}


def convert_to_ruediger_dict(d):

    new_d = dict()

    for key in d:
        if key in ruediger_conversion_table:
            new_key = ruediger_conversion_table[key]
            new_d[new_key] = d[key]
        else:
            print(f"Unconverted key: {key}")

    return new_d


def convert_from_ruediger_dict(d):

    inv_map = {v: k for k, v in ruediger_conversion_table.items()}

    new_d = dict()

    for key in d:
        if key in inv_map:
            new_key = inv_map[key]
            new_d[new_key] = d[key]
        else:
            print(f"Unconverted key: {key}")

    return new_d


if __name__ == "__main__":
    pass
