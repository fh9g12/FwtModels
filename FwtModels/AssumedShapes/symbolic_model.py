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
from .numeric_model import NumericModel


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

    def __init__(self,generalisedCoords,z_w,alpha_w,z_t,alpha_t,FwtParams,
            thetaIndex = None):
        self.p = FwtParams
        self.thetaIndex = thetaIndex
        # symbols for aero forces       
        self._w_g, self._alpha = me.dynamicsymbols('w_g alpha')
        self._extraTerms = sym.Matrix([self.p.alpha_0,self._w_g])

        # create generalised coordinates
        me.mechanics_printing()

        self.q = generalisedCoords
        self.qd = self.q.diff(t)
        self.qdd = self.qd.diff(t)

        self.z_w = z_w
        self.alpha_w = alpha_w
        self.z_t = z_t
        self.alpha_t = alpha_t

        # shape functions on the flexural axis (for aero forces)
        self.kappa_w = self.z_w.subs(x,self.p.x_f)
        self.kappa_t = self.z_t.subs(x,self.p.x_f)

        # create lambda equation for each z/alpha component

    @classmethod
    def B1_T0_RLFwt(cls,FwtParams):
        """ create an instance of a Fwt model with:
            B1 - 1 Bending Shape
            T0 - 0 torsional Shapes
            RLFwt - Rigid Locked FWT 
        """
        p = FwtParams
        q_0 = me.dynamicsymbols('q_1')
        q = sym.Matrix([q_0])

        z_w = q_0*y**2
        alpha_w = sym.Rational(0,1)

        z_t = z_w.subs(y,p.s_w) + z_w.diff(y).subs(y,p.s_w)*y
        alpha_t = sym.Rational(0,1)
        return cls(q,z_w,alpha_w,z_t,alpha_t,p)

    def Z(self,qs,xVal,yVals):
        # sub in all the variables to z and alpha eqns
        params = dict((s,s.value) for s in self.p.GetTuple())
        z_wEq = sym.lambdify((self.q,x,y),self.z_w.subs(params))
        z_tEq = sym.lambdify((self.q,x,y),self.z_t.subs(params))

        return np.where(yVals<=self.p.s_w.value,
                z_wEq(qs,xVal,yVals),
                z_tEq(qs,xVal,yVals))
    
    def Zpd(self,qs,xVal,yDf,name = 'Z',tipPositive = True):
        # sub in all the variables to z and alpha eqns
        ys = yDf.index.to_numpy()
        params = dict((s,s.value) for s in self.p.GetTuple())
        z_wEq = sym.lambdify((self.q,x,y),self.z_w.subs(params))
        z_tEq = sym.lambdify((self.q,x,y),self.z_t.subs(params))

        z = np.where(ys<=self.p.s_w.value,
                z_wEq(qs,xVal,ys),
                z_tEq(qs,xVal,ys))
        if z[-1]<0 & tipPositive:
            z = z*-1
        yDf[name] = z
        return yDf

    def Alpha(self,qs,xVal,yVals):
        # sub in all the variables to z and alpha eqns
        params = dict((s,s.value) for s in self.p.GetTuple())
        alpha_wEq = sym.lambdify((self.q,x,y),self.alpha_w.subs(params))
        alpha_tEq = sym.lambdify((self.q,x,y),self.alpha_t.subs(params))

        return np.where(y<=self.p.s_w.value,
                alpha_wEq(qs,xVal,yVals),
                alpha_tEq(qs,xVal,yVals))


    def GenerateEoM(self):
        """ Generates the EoM for the system"""
        # Get Lagranagian
        U = self.GeneratePotentialEnergy()
        T = self.GenerateKineticEnergy()
        L = sym.Matrix([T-U])

        # solve lagrangian (LHS)
        term_1 = L.jacobian(self.qd).diff(t).T
        term_2 = L.jacobian(self.q).T
        LHS = term_1-term_2

        # Extract the mass matrix (using the fact highest order term is 
        # second derivative)
        self.M = sym.simplify(LHS.jacobian(self.qdd))

        # use the Mass Matrix to find the remainder of the LHS 
        self.K = sym.simplify(LHS - self.M*self.qdd).jacobian(self.q)

        self.B_w, self.C_w, self.G_w = self.GetGeneralisedWingForces()

        self.B_t, self.C_t, self.G_t = self.GetGeneralisedTipForces()

        self.B = sym.simplify(self.B_w + self.B_t)
        self.C = sym.simplify(self.C_w + self.C_t)
        self.G = sym.simplify(self.G_w + self.G_t)



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

    def DecomposeQ(self,Q):
        # get the B, C & G matrices
        remainder = Q
        B = Q.jacobian(self.qd)
        remainder = (Q - B*self.qd)
        C = remainder.jacobian(self.q)
        remainder = remainder - C*self.q
        G = sym.simplify((Q-B*self.qd) - C*self.q).jacobian(self._extraTerms)
        return (B,C,G)


    def GeneratePotentialEnergy(self):
        """ Returns the symbolic expression the represents the potential energy
        of the system """
        # potential energy stored in main wing from bend and twisting
        U = sym.Rational(1,2)*((self.z_w.diff(y,y)**2*self.p.EI)
            .integrate((y,0,self.p.s_w)).integrate((x,0,self.p.c)))
        U = U + sym.Rational(1,2)*((self.alpha_w.diff(y)**2*self.p.GJ)
            .integrate((y,0,self.p.s_w)).integrate((x,0,self.p.c)))

        # potential energy stored in hinge spring ( assume last generalised
        # coord in theta)
        if self.thetaIndex is not None:
            U = U + sym.Rational(1,2)*self.p.k_theta*self.q[-1]**2

        # potential energy stored in main wing from gravitational forces
        U = U + ((self.z_w*self.p.g*self.p.m_w)
            .integrate((x,0,self.p.c),(y,0,self.p.s_w)))
        U = U + ((self.z_t*self.p.g*self.p.m_t)
            .integrate((x,0,self.p.c),(y,0,self.p.s_t)))
        return U

    def GenerateKineticEnergy(self):
        """ Returns the symbolic expression the represents the kinetic energy 
        of the system """
        T = ((self.z_w.diff(t)**2*sym.Rational(1,2)*self.p.m_w)
            .integrate((x,0,self.p.c),(y,0,self.p.s_w)))
        T = T + ((self.z_t.diff(t)**2*sym.Rational(1,2)*self.p.m_t)
            .integrate((x,0,self.p.c),(y,0,self.p.s_t)))
        return T

    def createNumericInstance(self, subs = None):
        """Method to create a simplified instance of the EoM:
        Inputs:
        Subs - a tuple of FwtVariables objects to sub vals for.
            If it is none (default) all FwtVariable in the local FwtParameters 
            instance will be substitued in the equations
            If it is provided only these variables will be substitued. 
            Allows parameter sweeps to be completed.
          """
        if subs is None:
            variables = self.p.GetTuple()
        else:
            variables = subs
        # v is passed to lambdify to use the version of the sqrt method that 
        # returns complex numbers
        v = [{'sqrt':np.lib.scimath.sqrt},'numpy']

        # create a function for each symbolic matrix stored in this instance
        M_eq = sym.lambdify(variables,self.M,v)
        K_eq = sym.lambdify(variables,self.K,v)
        B_eq = sym.lambdify(variables,sym.simplify(self.B_w+self.B_t),v)
        C_eq = sym.lambdify(variables,sym.simplify(self.C_w+self.C_t),v)
        G_eq = sym.lambdify(variables,sym.simplify(self.G_w+self.G_t),v)

        params = dict((s,s.value) for s in self.p.GetTuple())
        z_wEq = sym.lambdify((self.q,x,y),self.z_w.subs(params))
        z_tEq = sym.lambdify((self.q,x,y),self.z_t.subs(params))
        alpha_wEq = sym.lambdify((self.q,x,y),self.alpha_w.subs(params))
        alpha_tEq = sym.lambdify((self.q,x,y),self.alpha_t.subs(params))

        Z_eq = lambda qs,xs,ys : np.where(ys<=self.p.s_w.value,z_wEq(qs,xs,ys),
                                        z_tEq(qs,xs,ys))
        Alpha_eq = lambda qs,xs,ys : np.where(ys<=self.p.s_w.value,
                                        alpha_wEq(qs,xs,ys),alpha_tEq(qs,xs,ys))

        # create a tuple of values to pass to the functions generated above
        vals = tuple(map(lambda x:x.value,variables))

        # return a NumericModel instance, passing in numeric matrices created 
        # from the functions generated above 
        return NumericModel(M_eq(*vals),K_eq(*vals),B_eq(*vals),C_eq(*vals),
                    G_eq(*vals),Z_eq,Alpha_eq,self.p)

