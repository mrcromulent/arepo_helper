from utilities import suppress_stdout_stderr
import scipy.optimize as opt
from const import msol
import ic


class WDUtils(object):

    core_comp   = {'xHe4': 0.0, 'xC12': 0.3, 'xO16': 0.7, 'xNe20': 0.0, 'xMg24': 0.0}
    temperature = 5e5  # Kelvin

    @classmethod
    def get_default_params(cls, params=None):
        keys = ["core_comp", "temperature"]
        d = dict()

        for k in keys:
            if params and k in params:
                d[k] = params[k]
            else:
                d[k] = getattr(cls, k)

        return d

    @classmethod
    def __get_mass_from_rho_c(cls, rho_c, eos, params, offset):
        return cls.get_mass_from_rho_c(rho_c, eos, params) - offset

    @classmethod
    def get_mass_from_rho_c(cls, rho_c, eos, params=None):

        d = cls.get_default_params(params)
        with suppress_stdout_stderr():
            wd = ic.create_wd(eos, rho_c, temp=d["temperature"], **d["core_comp"])

        return wd['dm'].sum() / msol

    @classmethod
    def get_rho_c_from_mass(cls, mass, eos, params=None):

        params = cls.get_default_params(params)
        multiplier = 10.0

        rho_guess   = 1e4  # g/cc
        mass_guess  = cls.get_mass_from_rho_c(rho_guess, eos, params)

        if not mass_guess > mass:
            while mass_guess <= mass:

                rho_guess   *= multiplier
                mass_guess  = cls.get_mass_from_rho_c(rho_guess, eos, params)

        a = 1/multiplier * rho_guess
        b = rho_guess
        print(f"Attempting to estimate rho_c in range [{'{:.5E}'.format(a)}, {'{:.5E}'.format(b)}] g/cc")

        # noinspection PyTypeChecker
        rho_c = opt.bisect(cls.__get_mass_from_rho_c, a, b, args=(eos, params, mass))

        return rho_c


if __name__ == "__main__":
    from pyhelm_eos import loadhelm_eos
    
    equation_of_state = loadhelm_eos("./eostable/helm_table.dat", "./eostable/species05.txt", True)
    central_density = WDUtils.get_rho_c_from_mass(0.55, equation_of_state)
    print(f"Central Density = {'{:.5E}'.format(central_density)}")
