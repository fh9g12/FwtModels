import sympy as sym
import numpy as np
import pandas as pd
from scipy.linalg import eig
import sympy.physics.mechanics as me
from sympy.physics.vector.printing import vpprint, vlatex
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
        f = sym.simplify(term_1-M*p.qdd) -term_2
        return cls(p,M,f,T,U,ExtForces)


    def __init__(self,FwtParams,M,f,T,U,ExtForces = None):
        """
        Initialise a Symbolic model of the form 
        $M\ddot{q}+f(\dot{q},q,t)-ExtForces(\dot{q},q,t) = 0$

        with the Symbolic Matricies M,f,and Extforces
        """
        p = FwtParams

        self.M = M
        self.f = f
        self.T = T
        self.U = U

        # set external force function
        if ExtForces == None:
            self.ExtForces = lambda tup,x,t:[0]*int(len(x)/2)
        else:
            self.ExtForces = ExtForces

        #generate lambda functions
        tup = p.GetTuple()
        # Mass Matrix Eqn
        self.M_func = sym.lambdify((tup,p.x),self.M,"numpy")
        #func eqn
        self.f_func = sym.lambdify((tup,p.x),self.f,"numpy")

        # potential energy function
        self.u_eqn = sym.lambdify((tup,p.x),self.U,"numpy")
        # kinetic energy function
        self.t_eqn = sym.lambdify((tup,p.x),self.T,"numpy")

    def subs(self,p,*args):
        """
        Creates a new instance of a Symbolic model with the substutions supplied
         in args applied to all the Matricies
        """
        ExtForces = self.ExtForces.subs(p,*args) if self.ExtForces is not None else None
        return SymbolicModel(p,self.M.subs(*args),self.f.subs(*args),
                            self.T.subs(*args),self.U.subs(*args),ExtForces)
     

    def deriv(self,t,x,tup):
        external = self.ExtForces(tup,x,t)
        accels = np.linalg.inv(self.M_func(tup,x))@(-self.f_func(tup,x)+external)

        state_vc = []
        for i in range(0,int(len(x)/2)):
            state_vc.append(x[(i)*2+1])
            state_vc.append(accels[i,0])
        return tuple(state_vc)

    def GetFixedPoint(self,t,p0,tup):
        from scipy.optimize import minimize
        v = [0]*len(p0)
        def ObjectiveFunc(p,t,tup):           
            x = [j for i in range(len(v)) for j in [p[i],v[i]]]
            res = self.deriv(t,x,tup)
            return sum((i**2 for i in res))
        res = minimize(ObjectiveFunc,p0,(t,tup))
        return res.x




    #calculate the total energy in the system
    def KineticEnergy(self,t,x,tup):
        return self.t_eqn(*tup,x)

    def PotentialEnergy(self,t,x,tup):
        return self.u_eqn(*tup,x)

    def Energy(self,t,x,tup):
        return self.KineticEnergy(t,x,tup) + self.PotentialEnergy(t,x,tup)