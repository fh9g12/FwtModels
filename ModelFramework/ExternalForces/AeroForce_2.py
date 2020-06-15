import sympy as sym
from . import ExternalForce
from ..helper_funcs import LineariseMatrix
import sympy.physics.mechanics as me

class AeroForce_2(ExternalForce):
    @classmethod
    def PerUnitSpan(cls,FwtParams,Transform,C_L,alphadot,M_thetadot,e,rootAlpha,alpha_zero = 0,include_drag=False):
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

        if include_drag:
            F_x = L_w*(sym.sin(ang) + sym.cos(ang)*p.ratio_DL)
            F_z = L_w*(sym.cos(ang) - sym.sin(ang)*p.ratio_DL)
        else:
            F_x = L_w*sym.sin(ang)
            F_z = L_w*sym.cos(ang)

        wrench = sym.Matrix([F_x,0,F_z,0,M_w,0])

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

    def msubs(self,*args):
        return AeroForce_2(me.msubs(self._Q,*args),me.msubs(self.dAlpha,*args))

    def integrate(self,*args):
        return AeroForce_2(self._Q.integrate(*args),self.dAlpha)




