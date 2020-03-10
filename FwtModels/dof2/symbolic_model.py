import sympy as sym
import numpy as np
import sympyTransforms as symt

import pandas as pd
from scipy.linalg import eig
import sympy.physics.mechanics as me
from dataclasses import dataclass, InitVar, field
from sympy.abc import x,y,t


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

        self.M = self.EoM.jacobian(p.qdd)
        self.f = sym.simplify(self.EoM-self.M*p.qdd)
        
        # get equations for each generalised coordinates acceleration
        self.a_eq = self.M**-1*-self.f

        # create func for each eqn (param + state then derivative as inputs)
        state_vc = []
        for i in range(0,p.qs):
            state_vc.append(p.qd[i])
            state_vc.append(self.a_eq[i,0])
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
            self.ExtForces = lambda p,x,t:[0]*p.qs
        else:
            self.ExtForces = ExtForces

    def deriv(self,t,x,tup):
        return tuple(i[0] for i in self.X_func(*tup,self.ExtForces(tup,x,t),x))
    
    #calculate the total energy in the system
    def KineticEnergy(self,x,FwtParams,t):
        tup = FwtParams.GetNumericTuple(x,t)
        return self.t_eqn(*tup,x)

    def PotentialEnergy(self,x,FwtParams,t):
        tup = FwtParams.GetNumericTuple(x,t)
        return self.u_eqn(*tup,x)

    def Energy(self,x,FwtParams,t):
        return self.KineticEnergy(x,FwtParams,t) + \
                self.PotentialEnergy(x,FwtParams,t)
    
    def CruiseAngleEqn(self,FwtParams, ExtForces = None):
        p =FwtParams
        # get forcing Matrix
        Q = self.ExtForces.Q() if ExtForces == None else ExtForces.Q()
        #make all velocities zero
        subs = {p.qd[0]:0,p.qd[1]:0,p.qdd[0]:0,p.qdd[1]:0}
        X_stationary = self.X.subs(subs)
        Q_stationary = Q.subs(subs)

        # sub Q into EoM
        EoM_stationary = X_stationary.subs({self.F[0]:Q_stationary[0],self.F[1]:Q_stationary[1]})

        # sub in all other values apart from V
        tup = p.GetTuple(ignore=['V'])
        eqs = sym.simplify(EoM_stationary.subs({v:v.value for v in tup}))[[1,3],:]

        # sub in a replacement for V**2
        V_sq = sym.Symbol('V_sq')
        eqs = eqs.subs({p.V**2:V_sq})

        # solve both equations for d so we can eliminate it
        y1 = sym.solve(eqs[0],p.q[1])[0]
        y2 = sym.solve(eqs[1],p.q[1])[0]

        # solve for velocity
        x = sym.solve(y1-y2,V_sq)[0]
        y = y1.subs(V_sq,x)

        
        # return the lambdified eqn
        return sym.lambdify((p.q[0]),[sym.sqrt(x),y])













