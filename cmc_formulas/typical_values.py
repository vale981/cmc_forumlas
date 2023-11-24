from . import ureg
from . import constants
from . import formulas
import types
import numpy as np


CPW = types.SimpleNamespace()
CPW.Z_r = 50 * ureg.ohm
CPW.v_0 = 1.3e8 * ureg.m / ureg.s
CPW.d = 1 * ureg.cm
CPW.l = formulas.c(CPW.v_0, CPW.Z_r)
CPW.c = formulas.l(CPW.v_0, CPW.Z_r)
CPW.C_r = formulas.C_r(CPW.v_0, CPW.Z_r, CPW.d)
CPW.L_r = formulas.L_r(CPW.v_0, CPW.Z_r, CPW.d)
CPW.ω_0 = formulas.ω_0(CPW.L_r, CPW.C_r)

Transmon = types.SimpleNamespace()
Transmon.ratio = 50
Transmon.E_C = 100 * ureg.megahertz * ureg.hbar * 2 * ureg.pi
Transmon.E_J = Transmon.ratio * Transmon.E_C
Transmon.C_Σ = formulas.C_Σ(Transmon.E_C)
Transmon.ω_q = formulas.ω_q(Transmon.E_C, Transmon.E_J)
Transmon.α = (Transmon.E_C / ureg.hbar).to(ureg.gigahertz)

V_Dmax = 100 * ureg.microvolt


class CoplanarWaveguide:
    def __init__(self, impedance=CPW.Z_r, speed_of_light=CPW.v_0, length=CPW.d):
        self._impedance = impedance
        self._speed_of_light = speed_of_light
        self._length = length

    @property
    def L(self):
        return formulas.L_r(self._speed_of_light, self._impedance, self._length)

    @property
    def C(self):
        return formulas.C_r(self._speed_of_light, self._impedance, self._length)

    @property
    def ω(self):
        return formulas.ω_0(self.L, self.C)

    @property
    def f(self):
        return formulas.ω_0(self.L, self.C) / (2 * np.pi)

    @property
    def Z(self):
        return self._impedance


class TransmonQubit:
    def __init__(self, E_C=Transmon.E_C, E_J=Transmon.E_J):
        self._charge_energy = E_C.to(ureg.gigahertz * ureg.hbar)
        self._junction_energy = E_J.to(ureg.gigahertz * ureg.hbar)

    def __repr__(self):
        return f"""
{self.__class__.__name__}({self.E_C}, {self.E_J})
        """

    def overview(self):
        return f"""ω:      {self.ω}
ω/2π:   {self.f}
E_J:    {self.E_J}
E_C:    {self.E_C}
f_J:    {self.f_J}
f_C:    {self.f_C}
ratio:  {self.ratio}"""

    @property
    def E_C(self):
        return self._charge_energy

    @property
    def f_C(self):
        return self._charge_energy.to(ureg.gigahertz * ureg.planck_constant)

    @E_C.setter
    def E_C(self, E_C):
        self._charge_energy = E_C.to(ureg.gigahertz * ureg.hbar)

    @property
    def E_J(self):
        return self._junction_energy

    @property
    def f_J(self):
        return self._junction_energy.to(ureg.gigahertz * ureg.planck_constant)

    @E_J.setter
    def E_J(self, E_J):
        self._junction_energy = E_J.to(ureg.gigahertz * ureg.hbar)

    @property
    def C_Σ(self):
        return formulas.C_Σ(self.E_C).to(ureg.femtofarad)

    @property
    def ratio(self):
        return self.E_J / self.E_C

    @property
    def ω(self):
        return formulas.ω_q(self._charge_energy, self._junction_energy)

    @property
    def f(self):
        return self.ω / (2 * np.pi)

    @classmethod
    def from_capacitances(cls, C_Σ=formulas.C_Σ(Transmon.E_C), E_J=Transmon.E_J):
        return cls(formulas.E_C(C_Σ), E_J)

    @classmethod
    def from_angular_frequency_and_charging_energy(
        cls, ω_q=Transmon.ω_q, E_C=Transmon.E_C
    ):
        return cls(E_C, formulas.E_J_from_ω(ω_q, E_C))

    @classmethod
    def from_angular_frequency_and_transmon_ratio(
        cls, ω_q=Transmon.ω_q, ratio=Transmon.ratio
    ):
        charge_energy = formulas.E_C_from_ω_and_ratio(ω_q, ratio)
        return cls(charge_energy, charge_energy * ratio)


def Δ(transmon: TransmonQubit, cpw: CoplanarWaveguide):
    return transmon.ω - cpw.ω


def g(transmon: TransmonQubit, cpw: CoplanarWaveguide, C_c):
    return formulas.g(cpw.ω, transmon.C_Σ, C_c, transmon.E_J, cpw.Z)


def g(transmon: TransmonQubit, cpw: CoplanarWaveguide, C_c):
    return formulas.g(cpw.ω, transmon.C_Σ, C_c, transmon.E_J, cpw.Z)


def dimless_mag(quantity):
    return quantity.to(ureg.dimensionless).magnitude


def dispersive_metrics(transmon: TransmonQubit, cpw: CoplanarWaveguide, C_c):
    return {
        "Anharmonicity, α/Δ": dimless_mag(
            transmon.E_C / abs(Δ(transmon, cpw) * ureg.hbar)
        ),
        "Coupling Strength, g/Δ": dimless_mag(
            g(transmon, cpw, C_c) / abs(Δ(transmon, cpw))
        ),
        "RWA counter-rotating": dimless_mag(
            abs(Δ(transmon, cpw)) / (transmon.ω + cpw.ω)
        ),
        "Coupling/Transmon Cap": dimless_mag(C_c / transmon.C_Σ),
        "Coupling/CPW Cap": dimless_mag(C_c / cpw.C),
    }
