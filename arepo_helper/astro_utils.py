from species import ArepoSpeciesList
from const import NA, KB
import numpy as np


alpha   = 0.66
beta    = 12.8
rho_min = 1e-11
rho_max = 2e14
e_min   = 1e11
e_max   = 9.3e20


def estimate_collinear_lagrange_point(m1, m2, d):
    return (m1 - np.sqrt(m1 * m2)) / (m1 - m2) * d


def p_wrt_e(rho, e, gamma=1.667):
    return e * rho * (gamma - 1)


def e_wrt_p(rho, p, gamma=1.667):
    return p / (rho * (gamma - 1))


def test_e(rho, e):
    if rho < rho_min or rho > rho_max or e < e_min or e > e_max:
        return False
    else:
        return np.log10(e) > alpha * np.log10(rho) + beta


def test_p(rho, p):
    return test_e(rho, e_wrt_p(rho, p))
