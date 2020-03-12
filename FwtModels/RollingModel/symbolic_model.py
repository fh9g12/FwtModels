import sympy as sym
import numpy as np
import sympyTransforms as symt

import pandas as pd
from scipy.linalg import eig
import sympy.physics.mechanics as me
from dataclasses import dataclass, InitVar, field
from sympy.abc import x,y,t

class SymbolicModel:

    def __init__(self,U,Elements,FwtParams, ExtForces = None):
        p = FwtParams
        # matrix for forcing terms        
        self.F = sym.Matrix(me.dynamicsymbols(f'f:{p.qs}'))

        # Calc K.E
        self.T = sym.Integer(1)
        # add K.E for each Rigid Element
        for ele in Elements:
            T_e = ele.CalcKE(p)
            self.T = self.T + T_e

        # calculate EoM
        self.U = U
        self.Lag = sym.Matrix([self.T-self.U])
        term_1 = self.Lag.jacobian(p.qd).diff(me.dynamicsymbols._t).T
        term_2 = self.Lag.jacobian(p.q).T
        self.EoM = sym.simplify(term_1-term_2)

        # add exteral forces to EoM
        self.EoM = self.EoM - self.F

        self.M = self.EoM.jacobian(p.qdd)
        self.f = sym.simplify(self.EoM-self.M*p.qdd)

        tup = p.GetTuple()
        # Mass Matrix Eqn
        self.M_func = sym.lambdify((*tup,p.x),self.M)
        #func eqn
        self.f_func = sym.lambdify((*tup,self.F,p.x),self.f)
        # potential energy function
        self.u_eqn = sym.lambdify((*tup,p.x),self.U)
        # kinetic energy function
        self.t_eqn = sym.lambdify((*tup,p.x),self.T)

        # set external force function
        if ExtForces == None:
            self.ExtForces = lambda p,x,t:[0]*p.qs
        else:
            self.ExtForces = ExtForces

    def deriv(self,t,x,tup):
        external = self.ExtForces(tup,x,t)
        accels = np.linalg.inv(self.M_func(*tup,x))@(-self.f_func(*tup,external,x))
        state_vc = []
        for i in range(0,int(len(x)/2)):
            state_vc.append(x[(i)*2+1])
            state_vc.append(accels[i,0])
        return tuple(state_vc)
    
    #calculate the total energy in the system
    def KineticEnergy(self,x,tup,t):
        return self.t_eqn(*tup,x)

    def PotentialEnergy(self,x,tup,t):
        return self.u_eqn(*tup,x)

    def Energy(self,x,tup,t):
        return self.KineticEnergy(x,tup,t) + \
                self.PotentialEnergy(x,tup,t)