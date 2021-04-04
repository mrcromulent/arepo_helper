from utilities import suppress_stdout_stderr
from scipy.interpolate import interp1d
from collections import OrderedDict
import utilities
from names import n
from createICs import create_particles_healpix
from ics import ArepoICs
from ic import create_wd
import scipy.interpolate as interpolate
import scipy.optimize as opt
from species import ArepoSpeciesList
from pyhelm_eos import loadhelm_eos
import matplotlib.pyplot as plt
from numpy import log10
from const import msol
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
        # fig.savefig("test_abundances.png", dpi=300)

    def plot_healpix(self, helm_table, species_file, **kwargs):
        data = self.convert_to_healpix(helm_table, species_file, **kwargs)

        x = data["pos"][:, 0]
        y = data["pos"][:, 1]
        z = data["pos"][:, 2]

        fig, ax = plt.subplots()
        ax.scatter(x, y, c=z, s=0.01)
        fig.show()
        fig.savefig("test.png", dpi=300)

    """
    The C function create_particles_healpix creates a 3D heapix distribution based on the 1D distribution given as input

    :param dict data: Dictionary containing the 1D distribution of the relevant variables. data should contain several 
    numpy arrays. 
     data["dm"]     - N x 1 float   - Changes in mass over shells
     data["r"]      - N x 1 float   - Radius of shell
     data["rho"]    - N x 1 float   - Density of shell
     data["count"]  - int           - N
     data["ncells"] - int           - N
     data["v"]      - N x 1 float   - Volume or velocity, not sure. You can make this one all zeros
     data["p"]      - N x 1 float   - Pressure
     data["csnd"]   - N x 1 float   - Speed of sound in shell
     data["u"]      - N x 1 float   - Internal energy
     data["xnuc"]   - 1 x nspecies float - Nuclear composition (also works for N x nspecies, but only considers the 
     first) row
    :param helm_eos_table eos: Object created by `loadhelm_eos`
    :param int npart: Number of particles (fine to leave as zero?)
    :param int nspecies: Number of species
    :param float pmass: Particle mass?
    :param bool makebox: 
    :param float boxsize: The boxsize (used to centre the particles)
    :param float boxres: Box resolution
    :param int randomizeshells: 
    :param int randomizeradii: 
    :param float minenergy:
    :param bool fixgridpressure:
    :param float griddensity:
    :param float boxfactor:
    :param float cx:
    :param float cy:
    :param float cz:
    :param float transition_pmass:
    :return: Dictionary with keys "pos", "mass", "rho", "u", "vel", "xnuc" and "count".
    :rtrype: dict
    """

    def convert_to_healpix(self, helm_table, species_file, boxsize, pmass,
                           makebox=True,
                           minenergy=1e14,
                           boxfactor=10.0,
                           boxres=32,
                           npart=0,
                           randomizeshells=True,
                           randomizeradii=False):

        wdmass = 0.55 * msol
        c12prop = 0.5
        o16prop = 1 - c12prop

        def wdCOgetMassFromRhoCExact(rhoc, eos, mass=0., temp=5e5):
            import ic
            wd = ic.create_wd(eos, rhoc, temp=temp, xC12=0.5, xO16=0.5)
            print(rhoc, wd['dm'].sum() / msol)
            return wd['dm'].sum() / msol - mass

        def wdCOgetRhoCFromMassExact(mass, eos, temp=5e5, xtol=1e2):
            rhoguess = 1e4
            massguess = wdCOgetMassFromRhoCExact(rhoguess, eos, temp) * msol

            if massguess > mass:
                print("Already a central density of 1e4 g/ccm produces a WD more massive than %g msun (%g)." % (
                mass / msol, massguess / msol))

            while massguess <= mass:
                rhoguess *= 10.
                massguess = wdCOgetMassFromRhoCExact(rhoguess, eos, temp=temp) * msol

            print("Guess for central density: %g g/ccm." % rhoguess)
            rhoc = opt.bisect(wdCOgetMassFromRhoCExact, 0.1 * rhoguess, rhoguess, args=(eos, mass / msol, temp,),
                              xtol=xtol)
            print("Actual central density for %g msun WD: %g g/ccm." % (mass / msol, rhoc))
            return rhoc



        eos = loadhelm_eos(helm_table, species_file, True)  # Unused (hopefully)
        rhoc = wdCOgetRhoCFromMassExact(wdmass, eos)

        wd = create_wd(eos, rhoc, xC12=c12prop, xO16=o16prop)
        frho = interpolate.interp1d(wd['r'], wd['rho'], kind='cubic')
        fpres = interpolate.interp1d(wd['r'], wd['p'], kind='cubic')

        sp = ArepoSpeciesList(species_file)
        wd['v'] = np.zeros(wd['ncells'])
        wd['xnuc'] = np.zeros(sp.num_species)
        wd['xnuc'][sp.index_of('C12')] = c12prop
        wd['xnuc'][sp.index_of('O16')] = o16prop
        wd['count'] = wd['ncells']

        # df = self.data
        # gamma = 5. / 3.
        #
        # wd = dict()
        #
        # # Fill in the missing middle
        #
        # # n_extra_points = 200
        # n_extra_points = 0
        # npoints = df.shape[0] + n_extra_points
        #
        # extra_r = np.linspace(0, df['r'][0], n_extra_points)
        # extra_p = np.linspace(df['P'][0], df['P'][0], n_extra_points)
        # extra_rho = np.linspace(df['Rho'][0], df['Rho'][0], n_extra_points)
        #
        # wd['count'] = npoints
        # wd['ncells'] = npoints
        # wd['v'] = np.zeros(npoints)
        # wd['r'] = np.insert(np.array(df['r']), 0, extra_r)
        # wd['p'] = np.insert(np.array(df['P']), 0, extra_p)
        # wd['rho'] = np.insert(np.array(df['Rho']), 0, extra_rho)
        #
        # wd['csnd'] = np.sqrt(gamma * wd['p'] / wd['rho'])
        # wd['u'] = wd['p'] / (wd['rho'] * (gamma - 1))
        #
        # dm = np.zeros(npoints, dtype=np.float64)
        # wdec_q = df['-log(1-Mr/M*)']
        # qframe = np.array(wdec_q)
        # qframe = np.delete(qframe, 0)
        # qframe = np.delete(qframe, 0)
        # qframe = np.insert(qframe, 0, np.linspace(0, wdec_q[1], n_extra_points + 1))

        # for index, q in enumerate(qframe):
        #
        #     if index > 0:
        #         q2 = q
        #         q1 = qframe[index - 1]
        #
        #         val2 = np.power(10.0, -q2, dtype=np.float64)
        #         val1 = np.power(10.0, -q1, dtype=np.float64)
        #         dm[index] = msol * (val1 - val2)
        #
        # # When computing dm, I found that it was short by approximately 0.1 solar masses so I have added this as a
        # # kludge to ensure the dms have approximately the right distribution. Essentially, I am scaling up the whole
        # # distribution to sum to the correct value.
        # wd['dm'] = self.gp["m_tot"] * msol * dm / dm.sum()
        #
        # wd['xnuc'] = np.zeros((1, sp.num_species))
        # wd['xnuc'][0, sp.index_of("He4")] = df["X_He"][0]
        # wd['xnuc'][0, sp.index_of("C12")] = df["X_C"][0]
        # wd['xnuc'][0, sp.index_of("O16")] = df["X_O"][0]
        #
        # print(f"Input radii: {wd['r']}")

        with suppress_stdout_stderr():
            data = create_particles_healpix(wd, eos, npart,
                                            minenergy=minenergy,
                                            boxfactor=boxfactor,
                                            boxres=boxres,
                                            randomizeshells=randomizeshells,
                                            randomizeradii=randomizeradii,
                                            makebox=makebox,
                                            nspecies=sp.num_species,
                                            boxsize=boxsize,
                                            pmass=pmass)

        # data['pos'] += 0.5 * boxsize
        # center = np.array([boxsize / 2, boxsize / 2, boxsize / 2])
        #
        # r = np.linalg.norm(data['pos'] - center, axis=1)
        # r.sort()
        # print(r)
        # allowable_wdec_names = ["H1", "He4", "C12", "O16"]

        # for species in sp.species_dict:
        #     print(f"Interpolating species: {species}")
        #     sp_index = sp.index_of(species)
        #
        #     if species not in allowable_wdec_names:
        #         f = lambda rad: 0.0
        #     else:
        #
        #         yy = np.array(df[self.convert_species_name(species)])
        #
        #         x = np.array(wd['r'])
        #         y = np.insert(yy, 0, np.linspace(yy[0], yy[0], n_extra_points))
        #         f = interp1d(x, y)
        #
        #         # x = df['r']
        #         # y = df[self.convert_species_name(species)]
        #         #
        #         # if species == "He4":
        #         #     final_val = 1.0
        #         # else:
        #         #     final_val = 0.0
        #         #
        #         # def f(rad):
        #         #     if rad >= r_star:
        #         #         return final_val
        #         #     else:
        #         #         return interp1d(x, y)(rad)
        #
        #     for i, radius in enumerate(r):
        #         data['xnuc'][i, sp_index] = f(radius)

        # data['pres'] = np.zeros(data['rho'].shape[0])
        # data['temp'] = np.zeros(data['rho'].shape[0])
        # nparticles = data['xnuc'].shape[0]
        # data["id"] = np.array(range(1, nparticles + 1))

        print(f"rmax = {wd['r'].max()}")

        rad = np.sqrt(((data['pos'] - 0.5 * boxsize) ** 2.).sum(axis=1))
        i, = np.where(rad < wd['r'].max())
        rho = frho(rad[i])
        pres = fpres(rad[i])
        xnuc = np.zeros(sp.num_species)
        xnuc[sp.index_of('C12')] = c12prop
        xnuc[sp.index_of('O16')] = o16prop

        # For the cells included in the WD, find the internal energy and temperature
        for index in range(np.size(i)):
            idx = i[index]
            temp, data['u'][idx] = eos.pgiven(rho[index], xnuc, pres[index])
            if (data['u'][idx] < 1):
                print(data['u'][idx])
            # data['pres'][idx] = pres[index]

        if not (np.all(data['pos'] > 0)):
            print("BREAK")

        # Set mass equal to density in the WD, zero everywhere else
        data['mass'][:] = 0.
        data['mass'][i] = rho

        return data

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

    def create_particles_cube(self):
        pass


if __name__ == '__main__':
    wdec_directory = "/home/pierre/wdec/"
    helm_directory = "./data/snapshots/helm_table.dat"
    spec_directory = "./data/eostable/species13.txt"
    boxsize = 1e10
    pmass = 1e-6 * msol

    test = WDECResults(wdec_directory)
    data = test.convert_to_healpix(helm_directory, spec_directory, boxsize, pmass)

    # abund_sum = np.apply_along_axis(np.sum, 1, data['xnuc'])
    # fig, ax = plt.subplots()
    # ax.plot(abund_sum)
    # fig.show()

    # nparticles = data['xnuc'].shape[0]
    #
    # particle_list = np.array([nparticles, 0, 0, 0, 0, 0], dtype=np.int32)
    # particle_dict = OrderedDict({"PartType0": utilities.convert_from_ruediger_dict(data)})
    # header = {n.NUMPARTTOTAL: particle_list,
    #           n.NUMPARTTHISFILE: particle_list,
    #           n.NUMPARTTOTALHIGHWORD: np.zeros(6, dtype=np.int32),
    #           n.FLAGDOUBLEPRECISION: np.array(True, dtype=np.int32),
    #           n.BOXSIZE: boxsize}

    # aic = ArepoICs("bin.dat.ic.hdf5")
    # aic.write_ics(particle_dict, header)
    #
    # print(f"Number of unique values: {len(set(data['xnuc'][:, 2]))}")

    def gadget_write_ics(filename, data, format='hdf5', transpose=True, time=0., skipxnuc=False, double=False,
                         longids=False, num_files=1, boxsize=0., verbose=False, masses=None):
        if format == 'gadget1':
            # gadget_write_ics_format1(filename, data, transpose, time, skipxnuc, double, longids, num_files, boxsize,
            #                          verbose, masses)
            pass
        elif format == 'gadget2':
            # gadget_write_ics_format2(filename, data, transpose, time, skipxnuc, double, longids, num_files, boxsize,
            #                          verbose, masses)
            pass
        elif format == 'hdf5':
            gadget_write_ics_format3(filename, data, time, double, longids, num_files, boxsize, masses)
        else:
            raise ValueError('Choose a valid file format')


    def gadget_write_ics_format3(filename, data, time, double, longids, num_files, boxsize, masses, skipxnuc=False):
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


    gadget_write_ics("bin.dat.ic", data, boxsize=boxsize)
