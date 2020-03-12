import sympy as sym
import numpy as np
import pandas as pd
from scipy.linalg import eig
import sympy.physics.mechanics as me
from sympy.physics.vector.printing import vpprint, vlatex
from sympy.utilities.codegen import codegen
from sympy.utilities.autowrap import autowrap
from dataclasses import dataclass, InitVar, field
from sympy.abc import x,y,t


class SymbolicModel:
    """
    An instance of a folding wing tip model using assumed shapes.

    Required inputs are:
        generalisedCoords - array of the generalised coordinate symbols 
            (must be dynamic symbols)
        z_w,z_t,alpha_w,alpha_t - sympy expressions of the z and alpha postion
            of the wing and FWT
        FwtParameters - instance of the FwtParameters class (with the symbols 
            used in the above expressions)
        thetaIndex - index of theta (hinge angle) in generalisedCoords 
            (so energy equation knows which one) if no theta coordinate leave
            as 'None'
    """

    def __init__(self,Elements,FwtParams, ExtForces = None):

        p = FwtParams
        # matrix for forcing terms        
        self.F = sym.Matrix(me.dynamicsymbols(f'f:{p.qs}'))

        # Calc K.E
        self.T = sym.Integer(1)
        # add K.E for each Rigid Element
        for ele in Elements:
            T_e = ele.CalcKE(p)
            self.T = self.T + T_e

        # calc P.E
        self.U = sym.Integer(1)
        # add K.E for each Rigid Element
        for ele in Elements:
            U_e = ele.CalcPE(p)
            self.U = self.U + U_e

        # calculate EoM
        self.Lag = sym.Matrix([self.T-self.U])
        term_1 = self.Lag.jacobian(p.qd).diff(me.dynamicsymbols._t).T.expand()
        term_2 = self.Lag.jacobian(p.q).T


        self.M = term_1.jacobian(p.qdd) # assuming all parts in term 1 contribute only to mass matrix
        self.f = sym.simplify(term_1-self.M*p.qdd) -term_2 - self.F

        #self.EoM = term_1-term_2

        # add exteral forces to EoM
        #self.EoM = self.EoM - self.F

        #self.M = self.EoM.jacobian(p.qdd)
        #
        # self.f = self.EoM-self.M*p.qdd

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
            self.ExtForces = lambda tup,x,t:[0]*int(len(x)/2)
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

    def GetGeneralisedWingForces(self):
        """Returns the B C And G matrices for the generalised forces acting 
        upon the main wing section"""

        # Calculate external forces on the main wing
        # Calc forces on main wing (-1/2 as lift acts in oppisite direction to
        # z axis)
        dL_w = sym.Rational(-1,2)*self.p.rho*self.p.V**2*self.p.c*self.p.a_w* \
            (self.alpha_w + self.kappa_w.diff(t)/self.p.V+self.p.alpha_0+
            self._w_g/self.p.V)
        dM_w = sym.Rational(1,2)*self.p.rho*self.p.V**2*self.p.c**2* \
            (self.p.e*self.p.a_w*(self.alpha_w + self.kappa_w.diff(t)/self.p.V+
            self.p.alpha_0)+self.p.Malphadot*self.alpha_w.diff(t)*
            self.p.c/(4*self.p.V)+self._w_g/self.p.V)

        Q = ((sym.Matrix([self.kappa_w])
            .jacobian(self.q).T*dL_w)
            .integrate((y,0,self.p.s_w)))
        Q = Q + ((sym.Matrix([self.alpha_w])
            .jacobian(self.q).T*dM_w)
            .integrate((y,0,self.p.s_w)))

        # get the B, C & G matrices
        return self.DecomposeQ(Q)

    def GetGeneralisedTipForces(self):
        """Returns the B C And G matrices for the generalised forces acting
        upon the FWT section"""

        # Calculate external forces on the main wing
        # Calc forces on main wing (-1/2 as lift acts in oppisite direction to
        # z axis)
        dL_w = sym.Rational(-1,2)*self.p.rho*self.p.V**2*self.p.c*self.p.a_t* \
            (self.alpha_t + self.kappa_t.diff(t)/self.p.V+self.p.alpha_0+
            self._w_g/self.p.V)
        dM_w = sym.Rational(1,2)*self.p.rho*self.p.V**2*self.p.c**2* \
            (self.p.e*self.p.a_t*
            (self.alpha_t + self.kappa_t.diff(t)/self.p.V+self.p.alpha_0)+
            self.p.Malphadot*self.alpha_t.diff(t)*self.p.c/(4*self.p.V)+
            self._w_g/self.p.V)

        Q = ((sym.Matrix([self.kappa_t])
            .jacobian(self.q).T*dL_w)
            .integrate((y,0,self.p.s_t)))
        Q = Q + ((sym.Matrix([self.alpha_t])
            .jacobian(self.q).T*dM_w)
            .integrate((y,0,self.p.s_t)))

        # get the B, C & G matrices
        return self.DecomposeQ(Q)

    def PrintEoM(self):
        ql = vlatex(self.q)
        qdl = vlatex(self.qd)
        qddl = vlatex(self.qdd)

        lhs = vlatex(self.M)+qddl+'+'+ vlatex(self.K)+qddl
        rhs = vlatex(self.B)+qdl+'+'+vlatex(self.C)+ql+'+'\
            +vlatex(self.G)+vlatex(self._extraTerms)
        return lhs+'='+rhs

