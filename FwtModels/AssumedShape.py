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


class FwtVariable(sym.Symbol):
    """ child class of the symbol class to pepper it with a numerical value"""
    def __init__(self,v,sStr):
        self.value = v
        super().__init__()
    def __new__(cls,v,sStr):
        return super().__new__(cls,sStr)
        


class FwtParameters:
    m_w: FwtVariable = FwtVariable(0,'m_w')
    m_t: FwtVariable = FwtVariable(0,'m_t')
    x_f: FwtVariable = FwtVariable(0,'x_f')
    s_w: FwtVariable = FwtVariable(0,'s_w')
    s_t: FwtVariable = FwtVariable(0,'s_t')
    c: FwtVariable = FwtVariable(0,'c')
    Lambda: FwtVariable = FwtVariable(0,'Lambda')
    EI: FwtVariable = FwtVariable(0,'EI')
    GJ: FwtVariable = FwtVariable(0,'GJ')
    k_theta : FwtVariable = FwtVariable(0,'k_theta')
    rho: FwtVariable = FwtVariable(0,'rho')
    V: FwtVariable = FwtVariable(0,'V')
    a_w : FwtVariable = FwtVariable(0,'a_w')
    a_t : FwtVariable = FwtVariable(0,'a_t')
    alpha_0 : FwtVariable = FwtVariable(0,'alpha_0')
    e : FwtVariable = FwtVariable(0,'e')
    Malphadot : FwtVariable = FwtVariable(0,'M_alphadot')
    g : FwtVariable = FwtVariable(0,'g')
        
    
    def __init__(self,m_w = 0,m_t = 0,x_f=0,s_w=0,s_t=0,c=0,Lambda=0,EI=0,GJ=0,
        k_theta=0,rho=0,V=0,a_w=0,a_t=0,alpha_0=0,e=0,Malphadot=0,g = 9.81):
        self.m_w.value = m_w
        self.m_t.value = m_t
        self.x_f.value = x_f
        self.s_w.value = s_w
        self.s_t.value = s_t
        self.c.value = c
        self.Lambda.value = Lambda
        self.EI.value = EI
        self.GJ.value = GJ
        self.k_theta.value = k_theta
        self.rho.value = rho
        self.V.value = V
        self.a_w.value = a_w
        self.a_t.value = a_t
        self.alpha_0.value = alpha_0
        self.e.value = e
        self.Malphadot.value = Malphadot
        self.g.value = g
    
    def GetTuple(self):
        return (self.m_w,self.m_t,self.x_f,self.s_w,self.s_t,self.c,
        self.Lambda,self.EI,self.GJ,self.k_theta,self.rho,self.V,
        self.a_w,self.a_t,self.alpha_0,self.e,self.Malphadot,self.g)
    
@dataclass
class NumericModel:
    M : np.array
    K : np.array
    B : np.array
    C : np.array
    G : np.array
    Z : any
    Alpha : any
    Params : FwtParameters

    def FreeVibrationVals(self):
        """Calculate the eigen vectors and values of the free vibration 
        problem
        """
        return eig(self.K,self.M)

    def FreeVibrationValsPd(self):
        """ takes the results of the free vibration problem and places them in
        a data frame"""
        eVals,eVecs = self.FreeVibrationVals()

        # convert vectors into a 2D array and append the Freq's on the front
        vecs = np.array(eVecs)
        data = np.insert(vecs,0,0,axis=1)
        data[:,0] = np.sqrt(np.real(eVals))

        # Make column names for the data frame 
        coords = self.M.shape[1]
        cols = ['Freq']
        for i in range(0,coords):
            cols.append(f'Comp {i}')
        # Create the Data Frame
        df = pd.DataFrame(data,columns=cols)

        # return a sorted list of modes with Nans removed
        return df.dropna().sort_values(by='Freq').reset_index(drop=True)

    def AeroVibrations(self):
        Mp = self.Mprime()
        Kp = self.Kprime()
        ev,eve = eig(Kp,Mp)
        return (ev,eve)

    def Mprime(self):
        coords = self.M.shape[1]
        Mprime = np.eye(coords*2)
        Mprime[coords:,coords:]=self.M
        return Mprime

    def Kprime(self):
        coords = self.M.shape[1]
        Kprime = np.zeros((coords*2,coords*2))
        Kprime[coords:,coords:]=self.B
        Kprime[0:coords,coords:]=np.eye(coords)
        Kprime[coords:,0:coords]=self.C-self.K
        return Kprime

    def MprimeSym(self):
        coords = self.M.shape[1]
        Mprime = sym.eye(coords*2)
        Mprime[coords:,coords:]=self.M
        return Mprime

    def KprimeSym(self):
        coords = self.M.shape[1]
        Kprime = sym.zeros(coords*2)
        Kprime[coords:,coords:]=self.B
        Kprime[0:coords,coords:]=np.eye(coords)
        Kprime[coords:,0:coords]=self.C-self.K
        return Kprime

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
        self.thetaIndex = None
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
        if subs is None:
            variables = self.p.GetTuple()
        else:
            variables = subs
        v = [{'sqrt':np.lib.scimath.sqrt},'numpy']
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

        vals = tuple(map(lambda x:x.value,variables))

        return NumericModel(M_eq(*vals),K_eq(*vals),B_eq(*vals),C_eq(*vals),
                    G_eq(*vals),Z_eq,Alpha_eq,self.p)

