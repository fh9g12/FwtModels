import sympy as sym
from . import ExternalForce
from ..LambdifyExtension import msub

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
        self.dAlpha_func = sym.lambdify((tup,p.x),dAlpha,"numpy")
        super().__init__(p,Q)       

    def GetAlpha(self,tup,x,t,**kwargs):
        return self.dAlpha_func(tup,x)

    def linearise(self,p):
        x_subs = {(p.x[i],p.fp[i]) for i in range(-1,-len(p.x)-1,-1)}
        Q = self.Q()
        Q_p = Q.subs(x_subs)
        dAlpha_p = self.dAlpha.subs(x_subs)
        for i,x in enumerate(p.x):
            Q_p += Q.diff(x).subs(x_subs)*(x-p.fp[i])
            dAlpha_p += self.dAlpha.diff(x).subs(x_subs)*(x-p.fp[i])
        return AeroForce(p,Q_p,dAlpha_p)




