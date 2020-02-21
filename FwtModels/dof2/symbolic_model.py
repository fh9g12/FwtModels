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

    def __init__(self,U,Elements,FwtParams,ExtForces = None):
        p = FwtParams

        # matrix for forcing terms
        if ExtForces is not None:
            self.ExtForces = ExtForces
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
        if ExtForces is not None:
            self.EoM = self.EoM - self.F
        
        # get equations for each generalised coordinates acceleration
        self.a_eq = sym.linsolve(list(self.EoM[:]),list(p.qdd))

        # create func for each eqn (param + state then derivative as inputs)
        self.funcs = []
        tup = p.GetTuple()
        for i in range(0,p.qs):
            if ExtForces is None:
                self.funcs.append(sym.lambdify((*tup,p.q,p.qd),self.a_eq.args[0][i]))
            else:
                self.funcs.append(sym.lambdify((*tup,self.F,p.q,p.qd),self.a_eq.args[0][i]))
        # calculate U And T eqns
        
        self.u_eqn = sym.lambdify((*tup,p.q,p.qd),self.U)

        self.t_eqn = sym.lambdify((*tup,p.q,p.qd),self.T)
    
    def deriv(self,t,y,FwtParams):
        p=FwtParams
        result = ()
        qs=[]
        qds=[]
        for i in range(0,p.qs):
            qs.append(y[i*2])
            qds.append(y[i*2+1])
        
        tup = FwtParams.GetNumericTuple()
        for i in range(0,p.qs):
            result = (result + (y[i*2+1],))
            if self.ExtForces is None:
                v = self.funcs[i](*tup,qs,qds)
            else:
                fr = self.ExtForces.Calc(FwtParams,qs,qds,t)
                v = self.funcs[i](*tup,fr,qs,qds)           
            result = (result + (v,))
        
        return result
    
    #calculate the total energy in the system
    def KineticEnergy(self,y,FwtParams):
        qs,qds = self._getStates(y)
        tup = FwtParams.GetNumericTuple()
        return self.t_eqn(*tup,qs,qds)

    def PotentialEnergy(self,y,FwtParams):
        qs,qds = self._getStates(y)
        tup = FwtParams.GetNumericTuple()
        return self.u_eqn(*tup,qs,qds)

    def Energy(self,y,FwtParams):
        return self.KineticEnergy(y,FwtParams) + self.PotentialEnergy(y,FwtParams)

    def _getStates(self,y):
        qs=[]
        qds=[]
        for i in range(0,round(len(y)/2)):
            qs.append(y[i*2])
            qds.append(y[i*2+1])
        return qs,qds

