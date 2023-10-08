import numpy as np
from . import ureg
from . import constants


def Z_r(c, l):
    """
    Typical impedance as function of capacitance ``c`` and inductance
    ``l`` (per unit length).
    """
    return np.sqrt(c / l).to(ureg.ohm)


def v_0(c, l):
    """The speed of light in a CPW as function of capacitance ``c``
    and inductance ``l`` (per unit length).
    """
    return (1 / np.sqrt(c * l)).to(ureg.meter / ureg.second)


def c(v_0, Z_r):
    """
    Capacitance per unit length as function of the speed of light
    ``v_0`` and typical impedance ``Z_r``..
    """
    return (1 / (Z_r * v_0)).to(ureg.picofarad / ureg.meter)


def l(v_0, Z_r):
    """
    Inductance per unit length as function of the speed of light
    ``v_0`` and typical impedance ``Z_r``..
    """
    return (Z_r / (v_0)).to(ureg.nanohenry / ureg.meter)


def C_r(v_0, Z_r, d):
    """
    Total capacitance of a CPW.  Like ``c`` but with an extra
    argument: the length ``d``.
    """
    return (c(v_0, Z_r) * d).to(ureg.picofarad)


def L_r(v_0, Z_r, d):
    """
    Total inductance of a CPW.  Like ``c`` but with an extra
    argument: the length ``d``.
    """
    return (l(v_0, Z_r) * d).to(ureg.nanohenry)


def ω_0(L_r, C_r):
    """
    Fundamental (angular) frequency of CPW as a function of the
    _total_ inductance ``L_r`` and capacitance ``C_r``.
    """
    return (np.pi / np.sqrt(L_r * C_r)).to(ureg.gigahertz)


def ω_q(E_C, E_J):
    """
    The qubit (angular) frequency as a function as capacitive energy
    ``E_C`` and junction energy ``E_J``.
    """
    return ((np.sqrt(8 * E_J * E_C) - E_C) / ureg.hbar).to(ureg.gigahertz)


def E_C(C_Σ):
    return (ureg.elementary_charge**2 / (2 * C_Σ)).to(ureg.microjoule)


def C_Σ(E_C):
    return (ureg.elementary_charge**2 / (2 * E_C)).to(ureg.picofarad)


# @ureg.wraps(
#     ureg.picofarad, (ureg.second, ureg.picofarad, ureg.microjoule, ureg.hertz, ureg.ohm)
# )
def C_c(T_12, C_Σ, E_J, ω_B, Z_R=50 * ureg.ohm):
    e_C = E_C(C_Σ)
    Δ = np.abs(ω_q(e_C, E_J) - ω_B)
    return (
        np.sqrt(
            (Δ * C_Σ**2 * ureg.R_k)
            / (4 * T_12 * ω_B**2 * Z_R)
            * np.sqrt(2 * e_C / E_J)
        )
    ).to(ureg.picofarad)


def Δ(E_C, E_J, ω_B):
    return np.abs(ω_q(E_C, E_J) - ω_B)


def J(C_c, C_Σ, E_J, ω_B, Z_R=50 * ureg.ohm):
    e_C = E_C(C_Σ)
    return (
        ((4 * np.pi * ω_B**2 * C_c**2 * Z_R) / (Δ * C_Σ**2 * ureg.R_k))
        * np.sqrt(E_J / (2 * e_C))
    ).to(ureg.gigahertz)


def g(ω_0, C_Σ, C_c, E_J, Z_R=50 * ureg.ohm):
    e_C = E_C(C_Σ)
    return (
        ω_0
        * C_c
        / C_Σ
        * (E_J / (2 * e_C)) ** (1 / 4)
        * np.sqrt(2 * np.pi * Z_R / ureg.R_K)
    ).to(ureg.gigahertz)


def C_d_max(E_C, E_J, V_max):
    return ((2 * E_C / E_J) ** (1 / 4) * ureg.elementary_charge / (2 * V_max)).to(
        ureg.picofarad
    )
