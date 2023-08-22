from . import ureg
from . import constants
from . import formulas
import types


CPW = types.SimpleNamespace()
CPW.Z_r = 50 * ureg.ohm
CPW.v_0 = 1.3e8 * ureg.m / ureg.s
CPW.d = 1 * ureg.cm
CPW.l = formulas.c(CPW.v_0, CPW.Z_r)
CPW.c =  formulas.l(CPW.v_0, CPW.Z_r)
CPW.C_r = formulas.C_r(CPW.v_0, CPW.Z_r, CPW.d)
CPW.L_r = formulas.L_r(CPW.v_0, CPW.Z_r, CPW.d)
CPW.ω_0 = formulas.ω_0(CPW.L_r, CPW.C_r)

Transmon = types.SimpleNamespace()
Transmon.E_C = 400 * ureg.megahertz * ureg.hbar
Transmon.E_J = 25 * ureg.gigahertz * ureg.hbar
Transmon.C_Σ = formulas.C_Σ(Transmon.E_C)
Transmon.ω_q = formulas.ω_q(Transmon.E_C, Transmon.E_J)
Transmon.α = (Transmon.E_C / ureg.hbar).to(ureg.gigahertz)

V_Dmax = 100 * ureg.microvolt
