from typing import Union
import numpy as np


alpha   = 0.66
beta    = 12.8
rho_min = 1e-11
rho_max = 2e14
e_min   = 1e11
e_max   = 9.3e20


def estimate_collinear_lagrange_point(m1: float,
                                      m2: float,
                                      d: float) -> float:
    """Estimates collinear Lagrange point

    :param m1: Mass 1
    :param m2: Mass 2
    :param d: Distance between masses

    :return: Lagrange point
    """
    return (m1 - np.sqrt(m1 * m2)) / (m1 - m2) * d


def p_wrt_e(rho: Union[float, np.ndarray],
            e: Union[float, np.ndarray],
            gamma: Union[float, np.ndarray] = 1.667) -> Union[float, np.ndarray]:
    """

    :param rho:
    :param e:
    :param gamma:
    :return:
    """
    return e * rho * (gamma - 1)


def e_wrt_p(rho: Union[float, np.ndarray],
            p: Union[float, np.ndarray],
            gamma: Union[float, np.ndarray] = 1.667) -> Union[float, np.ndarray]:
    """

    :param rho:
    :param p:
    :param gamma:
    :return:
    """
    return p / (rho * (gamma - 1))


def test_e(rho: float,
           e: float) -> bool:
    """

    :param rho:
    :param e:
    :return:
    """
    if rho < rho_min or rho > rho_max or e < e_min or e > e_max:
        return False
    else:
        return np.log10(e) > alpha * np.log10(rho) + beta


def test_p(rho: float,
           p: float) -> bool:
    """

    :param rho:
    :param p:
    :return:
    """
    return test_e(rho, e_wrt_p(rho, p))
