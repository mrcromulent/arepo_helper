import matplotlib.pyplot as plt
import numpy as np


class ArepoEnergyFile:
    """ class for reading in energy.txt files """

    ntypes_dim = ["eins", "epots", "ekins", "etots", "masses"]

    def __init__(self, snappath="./"):

        if not snappath.endswith("/"):
            snappath += "/"
        self.snappath = snappath
        self.read_data()
        self.clean_data()

    def read_data(self):
        f = open(self.snappath + "energy.txt", "r")
        e = np.loadtxt(f)
        f.close()



        self.time = e[:, 0]

        self.length = len(self.time)

        self.ein = e[:, 1]
        self.epot = e[:, 2]
        self.ekin = e[:, 3]
        self.etot = e[:, 1:4].sum(axis=1)
        self.eins = np.zeros((6, len(self.time)))
        self.epots = np.zeros((6, len(self.time)))
        self.ekins = np.zeros((6, len(self.time)))
        self.etots = np.zeros((6, len(self.time)))
        self.masses = np.zeros((6, len(self.time)))
        for i in range(1, 7):
            self.eins[i - 1] = e[:, 1 + 3 * i]
            self.epots[i - 1] = e[:, 2 + 3 * i]
            self.ekins[i - 1] = e[:, 3 + 3 * i]
            self.etots[i - 1] = e[:, 1 + 3 * i:4 + 3 * i].sum(axis=1)
            self.masses[i - 1] = e[:, 21 + i]
        self.eadd = e[:, e.shape[1] - 1]

    def clean_data(self):
        """ remove overlapping intervals in time """
        import copy
        enew = copy.deepcopy(self)
        i = 1
        ilast = 0
        jlast = 0
        while i < len(self.time):
            if self.time[i] <= self.time[i - 1]:
                for j in range(i + 1, len(self.time)):
                    if self.time[j] > self.time[i - 1]:
                        # print("Time shift between %f d and %f d!"%(self.time[i]/86400., self.time[j]/86400.))
                        for k in self.__dict__.keys():
                            if not isinstance(self.__dict__[k], np.ndarray):
                                continue
                            if ilast == 0:
                                if len(enew.__dict__[k].shape) == 1:
                                    enew.__dict__[k] = enew.__dict__[k][:i]
                                else:
                                    tmp = []
                                    for m in range(self.__dict__[k].shape[0]):
                                        tmp.append(enew.__dict__[k][m][:i])
                                    enew.__dict__[k] = np.array(tmp)
                            else:
                                if len(enew.__dict__[k].shape) == 1:
                                    enew.__dict__[k] = np.concatenate((enew.__dict__[k], self.__dict__[k][jlast:i]))
                                else:
                                    tmp = []
                                    for m in range(self.__dict__[k].shape[0]):
                                        tmp.append(np.concatenate((enew.__dict__[k][m], self.__dict__[k][m][jlast:i])))
                                    enew.__dict__[k] = np.array(tmp)
                        ilast = i
                        jlast = j
                        break
                i = j
            else:
                i += 1
        if jlast > 0:
            for k in self.__dict__.keys():
                if not isinstance(self.__dict__[k], np.ndarray):
                    continue
                if len(enew.__dict__[k].shape) == 1:
                    enew.__dict__[k] = np.concatenate((enew.__dict__[k], self.__dict__[k][jlast:]))
                else:
                    tmp = []
                    for m in range(self.__dict__[k].shape[0]):
                        tmp.append(np.concatenate((enew.__dict__[k][m], self.__dict__[k][m][jlast:])))
                    enew.__dict__[k] = np.array(tmp)
        for k in self.__dict__.keys():
            self.__dict__[k] = enew.__dict__[k]

    def plot(self):
        fig, ax = plt.subplots()

        for k in self.__dict__.keys():
            if k not in ['snappath', 'time', 'length', 'ntypes_dim']:
                if k in self.ntypes_dim:
                    ax.plot(self.time, getattr(self, k)[0, :], label=f"{k}")
                else:
                    ax.plot(self.time, getattr(self, k), label=f"{k}")

        plt.legend()
        fig.show()


if __name__ == '__main__':
    test = ArepoEnergyFile("/run/user/1000/gvfs/sftp:host=gadi.nci.org.au/scratch/y08/ub0692/IsoWDs/wd0_35He/output/")
    test.plot()
