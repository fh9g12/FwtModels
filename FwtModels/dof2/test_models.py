from . import FwtParameters,FwtVariable
import sympy as sym
import numpy as np

def Simplified2DofModel(inner_freq,inner_mass,mass_ratio):
    # create the variables
    p = FwtParameters.Default2DoF() # parameter container

    p.m_w = FwtVariable(inner_mass,'m_w')   # mass of inner wing

    # set values for independent symbols
    p.s.value = 1
    p.c.value = 0.15
    p.g.value = 9.81
    p.rho.value = 1.225
    p.a_t.value = 2 * np.pi
    p.alpha_r.value = np.deg2rad(5)
    p.Lambda.value = np.deg2rad(10)

    # set expression for dependent symbols
    p.k = (sym.Integer(inner_freq)*2*sym.pi)**2*(p.m_w+p.m_w/sym.Float(mass_ratio))
    p.m = p.m_w/sym.Integer(mass_ratio)
    p.I_xx = sym.Rational(1,12)*p.m*p.s**2*1
    p.l = p.s*sym.Rational(1,2)
    return p