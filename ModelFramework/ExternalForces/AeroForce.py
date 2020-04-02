import sympy as sym
from . import ExternalForce
from ..LambdifyExtension import msub
from ..helper_funcs import LineariseMatrix

class AeroForce(ExternalForce):

    @classmethod
    def Generate(cls,FwtParams,Transform,C_L,int_tuple,alphadot,M_thetadot,e,rootAlpha,deltaAlpha,alpha_zero = 0):
        p = FwtParams
        ## force per unit length will following theredosons pseado-steady theory

        # add z velocity due to motion
        v_z_eff = sym.simplify(Transform.BodyVelocity()[2])

        # combine to get effective AoA
        dAlpha = alpha_zero + rootAlpha + deltaAlpha - v_z_eff/p.V

        # Calculate the lift force
        dynamicPressure = sym.Rational(1,2)*p.rho*p.V**2
        dL_w = dynamicPressure*p.c*C_L*dAlpha

        # Calulate the pitching Moment
        dM_w = dL_w*e*p.c # Moment due to lift
        dM_w += dynamicPressure*p.c**2*(M_thetadot*alphadot*p.c/(sym.Integer(4)*p.V))

        ## joint torques for lift are calculated in a frame aligned with the chordwise velocity direction
        wrench_lift = sym.Matrix([0,0,dL_w,0,0,0])
        velocity_frame = Transform

        dQ_L = (velocity_frame.ManipJacobian(p.q)).T
        dQ_L *= velocity_frame.InvAdjoint().T
        dQ_L = sym.trigsimp(dQ_L) # generally some sin^2 + cos^2 in here
        dQ_L *= wrench_lift

        ## joint torques for lift are calculated in a frame aligned with the chordwise velocity direction
        wrench_moment = sym.Matrix([0,0,0,0,dM_w,0])
        velocity_frame = Transform.R_y(deltaAlpha)

        dQ_M = (velocity_frame.ManipJacobian(p.q)).T
        dQ_M *= velocity_frame.InvAdjoint().T
        dQ_M *= wrench_moment

        dQ = dQ_L + dQ_M

        # generalised force
        _Q = dQ.integrate(int_tuple)
        return cls(p,_Q,dAlpha)
        
    def __init__(self,p,Q,dAlpha):
        tup = p.GetTuple()
        self.dAlpha = dAlpha
        self.dAlpha_func = sym.lambdify((tup,p.x,p.y_1),dAlpha,"numpy")
        super().__init__(p,Q)       

    def GetAlpha(self,tup,x,t,y):
        return self.dAlpha_func(tup,x,y)

    def linearise(self,p):
        Q_lin = LineariseMatrix(self.Q(),p.x,p.fp)
        dAlpha_lin = LineariseMatrix(self.dAlpha,p.x,p.fp)
        return AeroForce(p,Q_lin,dAlpha_lin)
    
    def subs(self,p,*args):
        return AeroForce(p,self._Q.subs(*args),self.dAlpha.subs(*args))




