import sympy as sym
from . import ExternalForce
from ..LambdifyExtension import msub
from ..helper_funcs import LineariseMatrix

class AeroForce_2(ExternalForce):

    @classmethod
    def PerUnitSpan(cls,FwtParams,Transform,C_L,alphadot,M_thetadot,e,rootAlpha,alpha_zero = 0):
        p = FwtParams
        ## force per unit length will following theredosons pseado-steady theory

        # add z velocity due to motion
        BodyJacobian = cls._trigsimp(Transform.BodyJacobian(p.q))

        v_z_eff = (BodyJacobian*p.qd)[2]
        
        # combine to get effective AoA
        dAlpha = alpha_zero + rootAlpha - v_z_eff/p.V

        # Calculate the lift force
        dynamicPressure = sym.Rational(1,2)*p.rho*p.V**2
        L_w = dynamicPressure*p.c*C_L*dAlpha

        # Calulate the pitching Moment
        M_w = L_w*e*p.c # Moment due to lift
        M_w += dynamicPressure*p.c**2*(M_thetadot*alphadot*p.c/(sym.Integer(4)*p.V))

        ## joint torques for lift are calculated in a frame aligned with the chordwise velocity direction
        ang = rootAlpha - v_z_eff/p.V
        wrench = sym.Matrix([L_w*sym.sin(ang),0,L_w*sym.cos(ang),0,M_w,0])

        _Q = BodyJacobian.T*wrench

        return cls(_Q,dAlpha)
        
    def __init__(self,Q,dAlpha):
        self.dAlpha = dAlpha
        super().__init__(Q) 

    @staticmethod
    def _trigsimp(expr):
        return sym.trigsimp(sym.powsimp(sym.cancel(sym.expand(expr))))      

    def linearise(self,x,x_f):
        Q_lin = LineariseMatrix(self.Q(),x,x_f)
        dAlpha_lin = LineariseMatrix(self.dAlpha,x,x_f)
        return AeroForce_2(Q_lin,dAlpha_lin)
    
    def subs(self,*args):
        return AeroForce_2(self._Q.subs(*args),self.dAlpha.subs(*args))

    def integrate(self,*args):
        return AeroForce_2(self._Q.integrate(*args),self.dAlpha)




