import matplotlib.pyplot as plt
from numpy import log10
import pandas as pd
import numpy as np
import os


class WDECResults(object):
    required_files = ["inputprof", "controlparams", "gridparameters", "corsico.dat", "output.dat"]
    names1 = ["n", "r", "Mr", "Lr", "T", "Rho", "P", "Eps", "Kap", "Cv"]
    names2 = ["n", "Chr", "Cht", "Epsr", "Epst", "Kapr", "Kapt", "Del", "Delad", "XHe"]
    names3 = ["n", "Fr", "log(q)", "ln(r/p)", "U", "V", "Voga1", "Ra", "chToR", "C1"]
    names4 = ["n", "C4", "B1", "B2", "B3", "B4", "alfa", "Tthl", "derdelad", "eta"]
    stpms = -4.365e-3

    def __init__(self, directory):

        if not os.path.exists(directory):
            raise OSError
        else:
            self.directory = os.path.abspath(directory)
            self.files = dict()
            self.data = None
            self.gp = None
            self.ip = None

            self.check_for_wdec_files()
            self.read_all_output_data()
            self.extract_input_params()

    def check_for_wdec_files(self):

        for filename in self.required_files:
            full_path = self.directory + "/" + filename
            assert os.path.exists(full_path)
            self.files[filename] = os.path.abspath(full_path)

    def extract_input_params(self):

        self.ip = self.read_inputprof()
        self.gp = self.read_gridparam()

    def read_abundances(self):
        fn = self.files["corsico.dat"]

        with open(fn, "r") as f:
            names = next(f)
            names = names.replace("#", "").split()
        df = pd.read_csv(fn, sep="\s+", names=names, skiprows=1)

        return df

    def read_gridparam(self):

        gridparam = dict()

        with open(self.files["gridparameters"], "r") as f:
            line = next(f)
            vals = list(line.split())
            num_vals = len(vals)

            if num_vals == 9:
                version = 15
            elif num_vals == 15:
                version = 16
            else:
                raise ValueError

            gridparam["version"] = version

            gridparam["t_eff"] = float(vals[0])
            gridparam["m_tot"] = float(vals[1]) / 1000
            gridparam["m_env"] = float(vals[2]) / 100
            gridparam["m_hel"] = float(vals[3]) / 100
            gridparam["m_hyd"] = float(vals[4]) / 100

            gridparam["he_frac"] = float(vals[5])

            gridparam["diff_coeff_hhe"] = float(vals[6])
            gridparam["diff_coeff_che"] = float(vals[7])

            if version == 15:
                gridparam["convective_efficiency"] = float(vals[8]) / 100
            elif version == 16:
                gridparam["alpha"] = float(vals[8])

                h1 = float(vals[9]) / 100
                h2 = round(float(vals[10]) / 100 * h1, 2)
                h3 = round(float(vals[11]) / 100 * h2, 2)
                gridparam["hvals"] = [h1, h2, h3]

                w1 = float(vals[12]) / 100
                w2 = float(vals[13]) / 100
                w3 = float(vals[14]) / 100
                gridparam["wvals"] = [w1, w2, w3]

        return gridparam

    def read_inputprof(self):
        inputprof = dict()

        with open(self.files["inputprof"], "r") as f:
            lines = f.readlines()
            inputprof["smoothing"] = self.fortran_float(lines[2])
            inputprof["buffer_inner"] = self.fortran_float(lines[4])
            inputprof["buffer_outer"] = self.fortran_float(lines[5])

        return inputprof

    def determine_pandas_df_keywords(self):
        fn = self.files["output.dat"]

        # Find number of points created by WDEC
        with open(fn, 'r') as f:
            next(f)
            next(f)
            line = next(f)
            npoints = int(line.replace('points', ''))

        pandas_kwargs = [{"sep": "\s+", "names": self.names1, "skiprows": 0, "nrows": npoints},
                         {"sep": "\s+", "names": self.names2, "skiprows": 0, "nrows": npoints},
                         # Further data frames have the centre and edge removed
                         {"sep": "\s+", "names": self.names3, "skiprows": 0, "nrows": npoints - 2},
                         {"sep": "\s+", "names": self.names4, "skiprows": 0, "nrows": npoints - 2}]

        # Determine the number of rows to skip by comparing with the column headers
        with open(fn, "r") as f:
            for line_no, line in enumerate(f):
                for i, names in enumerate([self.names1, self.names2, self.names3, self.names4]):

                    # The column header 'log q' caused an extra whitespace split, so it's been replaced with log(q)
                    tmp = line.replace("log q", "log(q)").split()
                    if tmp == names:
                        pandas_kwargs[i]["skiprows"] = line_no + 1

        return pandas_kwargs

    def read_all_output_data(self):
        fn = self.files["output.dat"]
        pandas_kwargs = self.determine_pandas_df_keywords()

        df1 = pd.read_csv(fn, **pandas_kwargs[0])
        df2 = pd.read_csv(fn, **pandas_kwargs[1])
        df3 = pd.read_csv(fn, **pandas_kwargs[2])
        df4 = pd.read_csv(fn, **pandas_kwargs[3])

        tmp1 = pd.merge(left=df1, right=df2, on='n')
        tmp2 = pd.merge(left=df3, right=df4, on='n')
        df = pd.merge(left=tmp1, right=tmp2, on='n')

        abundances = self.read_abundances()
        self.data = pd.concat([df, abundances], axis=1)

    def find_region_12_boundary(self):
        menv_mr = self.mr(self.gp['m_env'])
        buffer_inner = self.ip["buffer_inner"]
        return self.q(menv_mr - buffer_inner)

    def find_region_23_boundary(self):
        m_h = self.gp['m_hyd']
        amhyhe_tmp = 10 ** (- m_h)
        amr_hyhe = self.mr(amhyhe_tmp)
        qhyhe = abs(log10(amr_hyhe))
        buffer_outer = self.ip["buffer_outer"]

        return qhyhe - buffer_outer

    def plot_results(self):

        df = self.data

        _, ax = plt.subplots(2, 2)
        df.plot(x="r", y=["P"], ax=ax[0, 0], logy=True, title="Pressure")
        df.plot(x="r", y=["T"], ax=ax[0, 1], logy=True, title="Temperature")
        df.plot(x="r", y=["Chr"], ax=ax[1, 0], title="Pressure gradient")
        df.plot(x="r", y=["Cht", "Del", "Delad"], ax=ax[1, 1], logy=True, title="Temperature gradients")
        plt.tight_layout()
        plt.show()

        _, ax = plt.subplots(2, 2)
        df.plot(x="r", y=["Rho"], ax=ax[0, 0], logy=True, title="Density")
        df.plot(x="r", y=["Kap"], ax=ax[0, 1], logy=True, title="Opacity")
        df.plot(x="r", y=["Cv"], ax=ax[1, 0], logy=True, title="Specific Heat Capacity")
        df.plot(x="r", y=["Kapr", "Kapt"], ax=ax[1, 1], title="Opacity gradients")
        plt.tight_layout()
        plt.show()

        _, ax = plt.subplots(2, 2)
        df.plot(x="r", y=["Eps", "Epsr", "Epst"], ax=ax[0, 0], title="Energy generation rate")
        df.plot(x="r", y=["Lr"], ax=ax[0, 1], title="Radiative luminosity")
        df.plot(x="r", y=["U", "V"], ax=ax[1, 0], logy=True)
        df.plot(x="r", y=["alfa"], ax=ax[1, 1])
        plt.tight_layout()
        plt.show()

    def plot_abundances(self):

        co = self.data
        h1, h2, h3 = self.gp['hvals']
        w1, w2, w3 = self.gp['wvals']
        region_12_boundary = self.find_region_12_boundary()
        region_23_boundary = self.find_region_23_boundary()

        fig, ax = plt.subplots()
        co.plot(x="-log(1-Mr/M*)", y=["X_O", "X_C", "X_He", "X_H"], ax=ax, color=['red', 'black', 'blue', 'green'])

        if self.gp["version"] == 16:
            x = np.array([0, w1, w1 + w2, w1 + w2 + w3, 10 ** self.stpms - 1e-2])
            y = np.array([h1, h1, h2, h3, 0])
            ax.plot(self.q(x), y, 'kx', label="HW points")

        ax.set_ylabel("Abundance")
        ax.set_xlabel("-log10(1-Mr/M*)")
        ax.set_xlim([0, 18])
        ax.axvline(x=region_12_boundary, label="RI/RII boundary", linestyle='--', color='k', alpha=0.6)
        ax.axvline(x=region_23_boundary, label="RII/RIII boundary", linestyle='--', color='k', alpha=0.6)
        ax.axvline(x=self.gp['m_env'], label="M_env", linestyle='--', color='y', alpha=0.3)
        ax.axvline(x=self.gp['m_hel'], label="M_he", linestyle='--', color='b', alpha=0.3)
        ax.axvline(x=self.gp['m_hyd'], label="M_h", linestyle='--', color='g', alpha=0.3)
        ax.legend()
        fig.show()

    @staticmethod
    def convert_species_name(species):
        """
        Converts ArepoSpeciesList indexes to WDEC indexes
        :param species:
        :return:
        """

        if species == "H1":
            return "X_H"
        elif species == "He4":
            return "X_He"
        elif species == "O16":
            return "X_O"
        elif species == "C12":
            return "X_C"
        else:
            raise ValueError

    @staticmethod
    def q(mr):
        return -log10(1 - mr)

    @staticmethod
    def mr(q):
        return 1.0 - 10 ** (-q)

    @staticmethod
    def mr_diff(q1, q2):
        """
        :param q1: Inner q value
        :param q2: Outer q value
        :return: The equivalent Delta Mr for the corresponding qs

        Delta Mr
        = Mr2 - Mr1
        = msol * (1 - 10 ** -q2 - (1 - 10 ** -q2))
        = msol *(10 ** -q1 - 10 ** -q2)
        But msol = 1.989 * 10 ** 33 in grams
        = 1.989 * (10 ** (-q1 + 33) - 10 ** (-q2 + 33))
        """
        return 1.989 * (10 ** (-q1 + 33) - 10 ** (-q2 + 33))

    @staticmethod
    def fortran_float(strnum):
        return float(strnum.replace("d-", "e-").replace("d", "e+"))

    @classmethod
    def check_values(cls, gridparam, inputprof):

        m_hyd = gridparam["m_hyd"]
        m_env = gridparam["m_env"]
        buffer_inner = inputprof["buffer_inner"]
        menv_mr = cls.mr(m_env)
        boundary = menv_mr - buffer_inner

        print(f"CHECK: menv_mr - buffer_inner > 0: {boundary > 0}")
        if not boundary > 0:
            raise ValueError

        print(f"CHECK: m_hyd > m_env + 2: {m_hyd >= m_env + 2}")
        if not m_hyd >= m_env + 2:
            raise ValueError

        if gridparam["version"] == 16:
            w1, w2, w3 = gridparam['wvals']
            print(f"CHECK: w1 + w2 + w3 < 0.95: {w1 + w2 + w3 <= 0.95}")
            if not w1 + w2 + w3 <= 0.95:
                raise ValueError

    def write_temp_profile(self):
        out_filename = "TemperatureProfil.txt"

        with open(out_filename, 'w') as f:

            n_part = len(self.data["r"])
            f.write(str(n_part - 3) + "\n")
            for i in range(n_part):
                r = self.data["r"][i]
                t = self.data["T"][i]
                f.write(f"{r} {t} \n")

    def __str__(self):

        gp = self.gp
        ip = self.ip

        h1, h2, h3 = gp['hvals']
        w1, w2, w3 = gp['wvals']

        s = \
            f"\nFILE: {self.files['inputprof']} \n" \
            f"Smoothing value: {ip['smoothing']}. \n" \
            f"buffer_inner: mr = {ip['buffer_inner']}. \n" \
            f"buffer_outer: q = {ip['buffer_outer']} \n" \
            f"\nFILE: {self.files['gridparameters']} \n" \
            f"01 - WD with effective temperature: {gp['t_eff']} [K]. \n" \
            f"02 - The total mass {gp['m_tot']} [MSol]. \n" \
            f"03 - The envelope mass is q = {gp['m_env']}. \n" \
            f"04 - The helium mass is q = {gp['m_hel']} \n" \
            f"05 - The hydrogen mass is q = {gp['m_hyd']} \n" \
            f"06 - The helium abundance in the C/He/H region is {gp['he_frac']} %. \n" \
            f"07 - The diffusion constant at the base of the envelope: {gp['diff_coeff_hhe']}. \n" \
            f"08 - The diffusion constant of He at the base of pure He: {gp['diff_coeff_che']}. \n" \
            f"09 - Convective Alpha: {gp['alpha']} \n" \
            f"10 - h1: {h1} \n" \
            f"11 - h2: {h2} \n" \
            f"12 - h3: {h3} \n" \
            f"13 - w1: {w1} \n" \
            f"14 - w2: {w2} \n" \
            f"15 - w3: {w3}"

        return s
