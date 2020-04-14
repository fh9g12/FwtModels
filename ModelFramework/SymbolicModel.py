import sympy as sym
import numpy as np
import pandas as pd
from scipy.linalg import eig
import sympy.physics.mechanics as me
from sympy.physics.vector.printing import vpprint, vlatex
from sympy.abc import x,y,t
from .LambdifyExtension import msub
from .helper_funcs import LineariseMatrix
from .NumericModel import NumericModel   

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
    @classmethod
    def FromElementsAndForces(cls,FwtParams,Elements,ExtForces = None):
        """
        Create a symbolic Model instance from a set Elements and external forces
        """
        p = FwtParams 

        # Calc K.E
        T = sym.Integer(0)
        # add K.E for each Rigid Element
        for ele in Elements:
            T_e = ele.CalcKE(p)
            T = T + T_e

        # calc P.E
        U = sym.Integer(0)
        # add K.E for each Rigid Element
        for ele in Elements:
            U_e = ele.CalcPE(p)
            U = U + U_e

        # calculate EoM
        Lag = sym.Matrix([T-U])
        term_1 = Lag.jacobian(p.qd).diff(me.dynamicsymbols._t).T.expand()
        term_2 = Lag.jacobian(p.q).T

        # Get Mass Matrix and 'internal' forcing term
        M = term_1.jacobian(p.qdd) # assuming all parts in term 1 contribute only to mass matrix
        f = sym.expand(term_1-M*p.qdd) -term_2
        return cls(M,f,T,U,ExtForces)


    def __init__(self,M,f,T,U,ExtForces = None):
        """Initialise a Symbolic model of the form 
        $M\ddot{q}+f(\dot{q},q,t)-ExtForces(\dot{q},q,t) = 0$

        with the Symbolic Matricies M,f,and Extforces
        """
        self.M = M
        self.f = f
        self.T = T
        self.U = U

        self.ExtForces = ExtForces

    def CreateNumericModel(self,p):
        return NumericModel(p,self.M,self.f,self.T,self.U,self.ExtForces)

    def subs(self,*args):
        """
        Creates a new instance of a Symbolic model with the substutions supplied
         in args applied to all the Matricies
        """
        ExtForces = self.ExtForces.subs(*args) if self.ExtForces is not None else None

        # handle zero kinetic + pot energies
        T = self.T if isinstance(self.T,int) else self.T.subs(*args)
        U = self.U if isinstance(self.U,int) else self.U.subs(*args)
        return SymbolicModel(self.M.subs(*args),self.f.subs(*args),
                            T,U,ExtForces)

    def linearise(self,p):
        """
        Creates a new instance of the symbolic model class in which the EoM have been 
        linearised about the fixed point p.q_0
        """
        # Calculate Matrices at the fixed point
        # (go in reverse order so velocitys are subbed in before positon)

        # get the full EoM's for free vibration and linearise
        eom = self.M*p.qdd + self.f
        eom_lin = LineariseMatrix(eom,p.x,p.fp)

        #extract linearised M
        M_lin = eom_lin.jacobian(p.qdd)

        #extract linerised f
        f_lin = (eom_lin - M_lin*p.qdd).doit().expand()

        # Linearise the External Forces
        extForce_lin = self.ExtForces.linearise(p.x,p.fp) if self.ExtForces is not None else None

        # create the linearised model and return it
        return SymbolicModel(M_lin,f_lin,0,0,extForce_lin)

    def ExtractMatrices(self,p):
        """
        From the current symbolic model extacts the classic matrices A,B,C,D,E as per the equation below
        A \ddot{q} + B\dot{q} + Cq = D\dot{q} + Eq
        """
        A = self.M
        B = self.f.jacobian(p.qd)
        C = self.f.jacobian(p.q)
        D = self.ExtForces.Q().jacobian(p.qd)
        E = self.ExtForces.Q().jacobian(p.q)
        return A,B,C,D,E

    def FreeBodyEigenProblem(self,p):
        """
        gets the genralised eigan matrices for the free body problem.
        They are of the form:
            |   I   0   |       |    0    I   |
        M=  |   0   M   |   ,K= |   -C   -B   |
        such that scipy.linalg.eig(K,M) solves the problem 

        THE SYSTEM MUST BE LINEARISED FOR THIS TO WORK
        """
        M = sym.eye(p.qs*2)
        M[-p.qs:,-p.qs:]=self.M

        K = sym.zeros(p.qs*2)
        K[:p.qs,-p.qs:] = sym.eye(p.qs)
        K[-p.qs:,:p.qs] = -self.f.jacobian(p.q)
        K[-p.qs:,-p.qs:] = -self.f.jacobian(p.qd)
        return K,M

    def GeneralEigenProblem(self,p):
        """
        gets the genralised eigan matrices for use in solving the frequencies / modes. 
        They are of the form:
            |   I   0   |       |    0    I   |
        M=  |   0   M   |   ,K= |   E-C  D-B  |
        such that scipy.linalg.eig(K,M) solves the problem 

        THE SYSTEM MUST BE LINEARISED FOR THIS TO WORK
        """
        M_prime = sym.eye(p.qs*2)
        M_prime[-p.qs:,-p.qs:]=self.M

        _Q = self.ExtForces.Q() if self.ExtForces is not None else sym.Matrix([0]*p.qs)

        f = (_Q-self.f)

        K_prime = sym.zeros(p.qs*2)
        K_prime[:p.qs,-p.qs:] = sym.eye(p.qs)
        K_prime[-p.qs:,:p.qs] = f.jacobian(p.q)
        K_prime[-p.qs:,-p.qs:] = f.jacobian(p.qd)

        return K_prime, M_prime