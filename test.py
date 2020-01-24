from FwtModels import AssumedShape
import numpy as np
import sympy as sym
import sympy.physics.mechanics as me
from sympy.abc import y,x,t
import pandas as pd

semiSpan = 1.345  # semi-span of the wing
ratioFwt = 0   # ratio of the wing semi span that is fwt
wingMass = 2.75   # mass of entire semi-span wing
chord = 0.15

mu = wingMass/(chord*semiSpan) # mass density of the wing

# with width b and hieght h
h = 0.005 # hieght of bar
b = 0.03 # width of bar
I_xx = b*h**3/12
I_yy = b**3*h/12
J = I_xx + I_yy

# Material properties
E = 190e9
G = 74e9

# create an instance of the class holding all the properties for the FWT
p = AssumedShape.FwtParameters(m_w = mu,
                        m_t = mu,
                        x_f = 0.25,
                        s_w = semiSpan*(1-ratioFwt),
                        s_t = semiSpan*ratioFwt,
                        c = chord,
                        Lambda = np.deg2rad(10),
                        EI = E*I_xx,
                        GJ = G*J,
                        k_theta = 0,
                        rho = 1.225,
                        V = 10,
                        a_w = 2*np.pi,
                        a_t = 2*np.pi,
                        alpha_0 = 0,
                        e = 0.25,
                        Malphadot = -1.2,
                        g = 9.81, 
                        ThetaLocked = True)
