import sympy as sym
import numpy as np
from scipy.linalg import eig
import sympy.physics.mechanics as me
from sympy.physics.vector.printing import vpprint, vlatex
from sympy.utilities.codegen import codegen
from sympy.utilities.autowrap import autowrap
from dataclasses import dataclass, InitVar, field


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
        
    
    def __init__(self,m_w = 0,m_t = 0,x_f=0,s_w=0,s_t=0,c=0,Lambda=0,EI=0,GJ=0,k_theta=0,rho=0,V=0,a_w=0,a_t=0,alpha_0=0,e=0,Malphadot=0,g = 9.81):
        self.m_w.v = m_w
        self.m_t.v = m_t
        self.x_f.v = x_f
        self.s_w.v = s_w
        self.s_t.v = s_t
        self.c.v = c
        self.Lambda.v = Lambda
        self.EI.v = EI
        self.GJ.v = GJ
        self.k_theta.v = k_theta
        self.rho.v = rho
        self.V.v = V
        self.a_w.v = a_w
        self.a_t.v = a_t
        self.alpha_0.v = alpha_0
        self.e.v = e
        self.Malphadot.v = Malphadot
        self.g.v = g
    
    def GetTuple(self):
        return (self.m_w,self.m_t,self.x_f,self.s_w,self.s_t,self.c,self.Lambda,self.EI,self.GJ,self.k_theta,
        self.rho,self.V,self.a_w,self.a_t,self.alpha_0,self.e,self.Malphadot,self.g)
    
@dataclass
class NumericModel:
    M : np.array
    K : np.array
    B : np.array
    C : np.array
    G : np.array
    Params : FwtParameters

    def FreeVibrationVals(self):
        return eig(self.K,self.M)

    def AeroVibrations(self):
        coords = self.M.shape[1]
        Mprime = np.eye(coords*2)
        Mprime[coords:,coords:]=self.M
        print(Mprime)

        Kprime = np.zeros((coords*2,coords*2))
        Kprime[coords:,coords:]=-self.B
        Kprime[0:coords,coords:]=np.eye(coords)
        Kprime[coords:,0:coords]=self.K-self.C
        print(Kprime)

        return eig(Kprime,Mprime)

class SymbolicModel:
    """
    an Instance of a folding wing tip model using assumed shapes.
    Requires the inputs:
    1. generalisedCoords - A column vector of the generalised coordinate (q)
    2. TransformationMatrix - a transformation matrix from q to (z_w,alpha_w,z_t,alpha_t)
    """

    def __init__(self,generalisedCoords):

        # properties of main wing
        self._m_w, self._m_t, self._x_f, self._s_w, self._s_t = sym.symbols('m_w m_t x_f s_w s_t')
        self._Lambda, self._c, self._EI, self._GJ, self._g, self._k_theta = sym.symbols('Lambda c EI GJ g k_theta')
        self._rho, self._V, self._a_w, self._a_t, self._alpha_0, self._e = sym.symbols('rho V a_w a_t alpha_0 e ')
        self._Malphadot = sym.symbols('M_alphadot')

        self._x, self._y, self._t = sym.symbols('x y t')
        # symbols for aero forces       
        self._w_g, self._alpha = me.dynamicsymbols('w_g alpha')

        self._extraTerms = sym.Matrix([self._alpha_0,self._w_g])

        # create generalised coordinates
        self._a, self._b, self._theta = me.dynamicsymbols('a b theta')
        me.mechanics_printing()

        self.q = sym.Matrix([self._a,self._b,self._theta])
        self.qd = self.q.diff(self._t)
        self.qdd = self.qd.diff(self._t)

        # Define shape functions
        # self.z_w = self._a*self._y**2+self._b*(self._x-self._x_f)*self._y
        # self.alpha_w = self._b*self._y
        self.z_w = self._a*self._y**2/self._s_w**2##+self._b*(self._x-self._x_f)*self._y
        self.alpha_w = sym.Rational(0,1)

        self.z_t = self.z_w.subs(self._y,self._s_w) + self._theta*(self._y)
        self.alpha_t = self.alpha_w.subs(self._y,self._s_w) + self._theta*sym.sin(self._Lambda)

        # shape functions on the flexural axis (for aero forces)
        self.kappa_w = self.z_w.subs(self._x,self._x_f)
        self.kappa_t = self.z_t.subs(self._x,self._x_f)

    def GenerateEoM(self):
        """ Generates the EoM for the system"""
        # Get Lagranagian
        U = self.GeneratePotentialEnergy()
        T = self.GenerateKineticEnergy()
        L = sym.Matrix([T-U])

        # solve lagrangian (LHS)
        term_1 = L.jacobian(self.qd).diff(self._t).T
        term_2 = L.jacobian(self.q).T
        LHS = term_1-term_2

        #Extract the mass matrix (using the fact highest order term is second derivative)
        self.M = sym.simplify(LHS.jacobian(self.qdd))

        # use the Mass Matrix to find the remainder of the LHS 
        self.K = sym.simplify(LHS - self.M*self.qdd).jacobian(self.q)

        self.B_w, self.C_w, self.G_w = self.GetGeneralisedWingForces()

        self.B_t, self.C_t, self.G_t = self.GetGeneralisedTipForces()

        self.B = sym.simplify(self.B_w + self.B_t)
        self.C = sym.simplify(self.C_w + self.C_t)
        self.G = sym.simplify(self.G_w + self.G_t)



    def GetGeneralisedWingForces(self):
        """Returns the B C And G matrices for the generalised forces acting upon the main wing section"""

        # Calculate external forces on the main wing
        # Calc forces on main wing (-1/2 as lift acts in oppisite direction to z axis)
        dL_w = sym.Rational(-1,2)*self._rho*self._V**2*self._c*self._a_w*(self.alpha_w + self.kappa_w.diff(self._t)/self._V+self._alpha_0+self._w_g/self._V)
        dM_w = sym.Rational(1,2)*self._rho*self._V**2*self._c**2*(self._e*self._a_w*(self.alpha_w + self.kappa_w.diff(self._t)/self._V+self._alpha_0)+self._Malphadot*self.alpha_w.diff(self._t)*self._c/(4*self._V)+self._w_g/self._V)

        Q = (sym.Matrix([self.kappa_w]).jacobian(self.q).T*dL_w).integrate((self._y,0,self._s_w))
        Q = Q + (sym.Matrix([self.alpha_w]).jacobian(self.q).T*dM_w).integrate((self._y,0,self._s_w))

        # get the B, C & G matrices
        return self.DecomposeQ(Q)

    def GetGeneralisedTipForces(self):
        """Returns the B C And G matrices for the generalised forces acting upon the FWT section"""

        # Calculate external forces on the main wing
        # Calc forces on main wing (-1/2 as lift acts in oppisite direction to z axis)
        dL_w = sym.Rational(-1,2)*self._rho*self._V**2*self._c*self._a_t*(self.alpha_t + self.kappa_t.diff(self._t)/self._V+self._alpha_0+self._w_g/self._V)
        dM_w = sym.Rational(1,2)*self._rho*self._V**2*self._c**2*(self._e*self._a_t*(self.alpha_t + self.kappa_t.diff(self._t)/self._V+self._alpha_0)+self._Malphadot*self.alpha_t.diff(self._t)*self._c/(4*self._V)+self._w_g/self._V)

        Q = (sym.Matrix([self.kappa_t]).jacobian(self.q).T*dL_w).integrate((self._y,0,self._s_t))
        Q = Q + (sym.Matrix([self.alpha_t]).jacobian(self.q).T*dM_w).integrate((self._y,0,self._s_t))

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
        """ Returns the symbolic expression the represents the potential energy of the system """
        # potential energy stored in main wing from bend and twisting
        U = sym.Rational(1,2)*(self.z_w.diff(self._y,self._y)**2*self._EI).integrate((self._y,0,self._s_w))
        U = U + sym.Rational(1,2)*(self.alpha_w.diff(self._y)**2*self._GJ).integrate((self._y,0,self._s_w))
        # potential energy stored in hinge spring
        U = U + sym.Rational(1,2)*self._k_theta*self._theta**2

        # potential energy stored in main wing from gravitational forces
        U = U + (self.z_w*self._g*self._m_w).integrate((self._x,0,self._c),(self._y,0,self._s_w))
        U = U + (self.z_t*self._g*self._m_t).integrate((self._x,0,self._c),(self._y,0,self._s_t))
        return U

    def GenerateKineticEnergy(self):
        """ Returns the symbolic expression the represents the kinetic energy of the system """
        T = (self.z_w.diff(self._t)**2*sym.Rational(1,2)*self._m_w).integrate((self._x,0,self._c),(self._y,0,self._s_w))
        T = T + (self.z_t.diff(self._t)**2*sym.Rational(1,2)*self._m_t).integrate((self._x,0,self._c),(self._y,0,self._s_t))
        return T

    def createNumericInstance(self,FwtParams : FwtParameters):
        variables = (self._m_w,self._m_t,self._x_f,self._s_w,self._s_t,self._c,self._Lambda,self._EI,self._GJ,self._k_theta,self._rho,self._V,self._a_w,self._a_t,self._alpha_0,self._e,self._Malphadot,self._g)
        M_eq = sym.lambdify(variables,self.M)
        K_eq = sym.lambdify(variables,self.K)
        B_eq = sym.lambdify(variables,sym.simplify(self.B_w+self.B_t))
        C_eq = sym.lambdify(variables,sym.simplify(self.C_w+self.C_t))
        G_eq = sym.lambdify(variables,sym.simplify(self.G_w+self.G_t))

        vals = FwtParams.GetTuple()

        return NumericModel(M_eq(*vals),K_eq(*vals),B_eq(*vals),C_eq(*vals),G_eq(*vals),FwtParams)

