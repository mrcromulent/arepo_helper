import matplotlib.pyplot as plt
from typing import Any, Dict
from .utilities import sci
import numpy as np
import warnings
import os


class ArepoEnergyFile:
    """Class to read AREPO energy.txt files."""

    fields = {
        "Time": [0],
        "Internal Energy": [1],
        "Potential Energy": [2],
        "Kinetic Energy": [3],
        "Internal Energies": [4, 9],
        "Potential Energies": [10, 15],
        "Kinetic Energies": [16, 21],
        "Total Masses": [22, 27],
        "Eadd": [28]
    }

    def __init__(self, directory):

        self.filepath = os.path.join(directory, "energy.txt")
        self.length = 0
        self.data = Dict[str, Any]
        self.read_data()
        self.sense_checks()

    def read_data(self):

        data = dict()

        with open(self.filepath, "r") as f:
            e = np.loadtxt(f)
            self.length = np.shape(e)[0]

        for field in self.fields:
            cols = self.fields[field]

            if len(cols) > 1:
                data[field] = e[:, cols[0]:cols[1]]
            else:
                data[field] = e[:, cols[0]]

        data["Total Energy"] = data["Internal Energy"] + \
                               data["Kinetic Energy"] - \
                               data["Potential Energy"]

        self.data = data

    def sense_checks(self):

        if np.any(self.data["Potential Energy"] == 0):
            warnings.warn("Potential Energy was not computed")

    def plot(self, savefig=False):
        fig, ax = plt.subplots()

        for field in self.data:
            this_data = self.data[field]
            one_d = len(np.shape(this_data)) == 1
            if one_d and field != "Time":
                ax.plot(self.data["Time"], this_data, label=f"{field}")

        max_val = max(self.data["Total Energy"])
        min_val = min(self.data["Total Energy"])
        max_pc_diff = abs((max_val - min_val) / min_val)

        ax.legend()
        ax.set_xlabel("Time [s]")
        ax.set_ylabel("Energy [erg]")
        ax.set_title(f"Energy balance. Max diff = {sci(max_pc_diff * 100)} %")

        if savefig:
            savepath = os.path.join(os.path.dirname(self.filepath), "energy_balance.png")
            fig.savefig(savepath, dpi=300)
        else:
            fig.show()
