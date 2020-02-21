import sympy as sym
import numpy as np
import sympyTransforms as symt

import pandas as pd
from scipy.linalg import eig
import sympy.physics.mechanics as me
from sympy.physics.vector.printing import vpprint, vlatex
from sympy.utilities.codegen import codegen
from sympy.utilities.autowrap import autowrap
from dataclasses import dataclass, InitVar, field
from sympy.abc import x,y,t
from .rigidElement import RigidElement, MassMatrix
import typing


class SymbolicModel:
    """
    An instance of a 2 Dof FWT model that is a pendulum attached to a 
    vertical spring.

    Required inputs are:
        generalisedCoords - array of the generalised coordinate symbols 
            (must be dynamic symbols)
        U - expression for the potenital energy
        Transform - transform for the single point mass
        FwtParameters - instance of the FwtParameters class (with the symbols 
            used in the above expressions)
    """

    def __init__(self,U,Elements,FwtParams, ExtForces = None):
        p = FwtParams

        # matrix for forcing terms
        # if ExtForces is not None:
        #     self.ExtForces = ExtForces
        
        self.F = sym.Matrix(me.dynamicsymbols(f'f:{p.qs}'))

        # Calc K.E
        self.T = sym.Integer(1)
        # add K.E for each Rigid Element
        for ele in Elements:
            T_e = ele.CalcKE(p.q,p.qd)
            self.T = self.T + T_e

        # calculate EoM
        self.U = U
        self.Lag = sym.Matrix([self.T-self.U])
        term_1 = self.Lag.jacobian(p.qd).diff(me.dynamicsymbols._t).T
        term_2 = self.Lag.jacobian(p.q).T
        self.EoM = sym.simplify(term_1-term_2)

        # add exteral forces to EoM
        self.EoM = self.EoM - self.F
        
        # get equations for each generalised coordinates acceleration
        self.a_eq = sym.linsolve(list(self.EoM[:]),list(p.qdd))

        # create func for each eqn (param + state then derivative as inputs)
        state_vc = []
        for i in range(0,p.qs):
            state_vc.append(p.qd[i])
            state_vc.append(self.a_eq.args[0][i])
        self.X = sym.Matrix(state_vc)
        
        tup = p.GetTuple()

        #state vector function
        self.X_func = sym.lambdify((*tup,self.F,p.x),self.X)
        # potential energy function
        self.u_eqn = sym.lambdify((*tup,p.x),self.U)
        # kinetic energy function
        self.t_eqn = sym.lambdify((*tup,p.x),self.T)

        # set external force function
        if ExtForces == None:
            self.__fr = lambda p,x,t:[0]*p.qs
        else:
            self.__fr = ExtForces

    
    def deriv(self,t,x,FwtParams):
        p=FwtParams
        tup = FwtParams.GetNumericTuple(x,t)
        return tuple(i[0] for i in self.X_func(*tup,self.__fr(p,x,t),x))
    
    #calculate the total energy in the system
    def KineticEnergy(self,x,FwtParams,t):
        tup = FwtParams.GetNumericTuple(x,t)
        return self.t_eqn(*tup,x)

    def PotentialEnergy(self,x,FwtParams,t):
        tup = FwtParams.GetNumericTuple(x,t)
        return self.u_eqn(*tup,x)

    def Energy(self,x,FwtParams,t):
        return self.KineticEnergy(x,FwtParams) + \
                self.PotentialEnergy(x,FwtParams)