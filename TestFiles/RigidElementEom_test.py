import sympy as sym
import sympy.physics.mechanics as me

import numpy as np
from scipy.integrate import odeint,LSODA,BDF,solve_ivp
import matplotlib.pyplot as plt

import sys, os

sys.path.insert(1, os.path.join(sys.path[0], '..'))
import sympyTransforms as symt
import custom_plot_objects as cpo
import FwtModels.dof2 as dof2

# create the variables
p = dof2.FwtParameters()
p.m.value = 3
p.l.value = 1
p.s.value = 2
p.c.value = 0.15
p.k.value = 100
p.g.value = 9.81

p.m_1 = p.AddParam(1.5,'m_1')
p.m_2 = p.AddParam(1.5,'m_2')


# get generalised coords (theta and d)
qs = 2
q = sym.Matrix(me.dynamicsymbols(f'q:{qs}'))

# potential Energy
half = sym.Rational(1,2)
U = p.m*p.g*(-p.l*sym.cos(q[0])) + half*p.k*(q[1])**2

# Create Point masses

# Transform for the mass


pendulum_frame = symt.HomogenousTransform().Translate(0,0,q[1]).R_x(-sym.pi*sym.Rational(1,2)+q[0])
tip_frame = pendulum_frame.Translate(0,p.l,0)
rot = sym.Matrix([q[0],0,0])

m1 = dof2.RigidElement.PointMass(pendulum_frame,rot,p.m_1)
m2 = dof2.RigidElement.PointMass(tip_frame,rot,p.m_2)

sm = dof2.SymbolicModel(q,U,[m1,m2],p)
