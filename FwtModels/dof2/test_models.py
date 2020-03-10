from .. import FwtParameters,FwtVariable
import sympy as sym
import numpy as np

def Simplified2DofModel(inner_freq,inner_mass,mass_ratio):
    # create the variables
    p = FwtParameters(2) # parameter container

    p.m: FwtVariable = FwtVariable(0,'m') # mass of FWT
    p.l: FwtVariable = FwtVariable(0,'l') # dist from hinge to CoM
    p.s: FwtVariable = FwtVariable(1,'s') # span
    p.c: FwtVariable = FwtVariable(0.15,'c') # chord
    p.k: FwtVariable = FwtVariable(0,'k') # spring constant
    p.g : FwtVariable = FwtVariable(9.81,'g') # gravity
    p.Lambda: FwtVariable = FwtVariable(np.deg2rad(10),'Lambda') # flare angle
    p.rho: FwtVariable = FwtVariable(1.225,'rho') # density
    p.V: FwtVariable = FwtVariable(0,'V') # velocity
    p.a_t : FwtVariable = FwtVariable(2 * np.pi,'a_t') # C_L of FWT
    p.alpha_r : FwtVariable = FwtVariable(np.deg2rad(5),'alpha_r') # C_L of FWT

    p.m_w = FwtVariable(inner_mass,'m_w')   # mass of inner wing

    # set expression for dependent symbols
    p.k = (sym.Integer(inner_freq)*2*sym.pi)**2*(p.m_w+p.m_w/sym.Float(mass_ratio))
    p.m = p.m_w/sym.Integer(mass_ratio)
    p.I_xx = sym.Rational(1,12)*p.m*p.s**2*1
    p.l = p.s*sym.Rational(1,2)
    return p