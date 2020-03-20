import sympy as sym
from . import ExternalForce

class AeroForce(ExternalForce):

    def __init__(self,FwtParams,Transform,C_L,int_tuple,rootAlpha,delta_alpha,alphadot,M_thetadot,e):
        p = FwtParams
        ## force per unit length will following theredosons pseado-steady theory

        # add z velocity due to motion
        v_z_eff = sym.simplify(Transform.BodyVelocity()[2])

        # combine to get effective AoA
        self.dAlpha = rootAlpha + delta_alpha - v_z_eff/p.V

        # Calculate the lift force
        dynamicPressure = sym.Rational(1,2)*p.rho*p.V**2
        self.dL_w = dynamicPressure*p.c*C_L*self.dAlpha

        # Calulate the pitching Moment
        self.dM_w = self.dL_w*e*p.c # Moment due to lift
        self.dM_w += dynamicPressure*p.c**2*(M_thetadot*alphadot*p.c/(sym.Integer(4)*p.V))

        ## joint torques for lift are calculated in a frame aligned with the chordwise velocity direction
        wrench_lift = sym.Matrix([0,0,self.dL_w,0,0,0])
        velocity_frame = Transform

        self.dQ_L = (velocity_frame.ManipJacobian(p.q)).T
        self.dQ_L *= velocity_frame.InvAdjoint().T
        self.dQ_L = sym.trigsimp(self.dQ_L) # generally some sin^2 + cos^2 in here
        self.dQ_L *= wrench_lift

        ## joint torques for lift are calculated in a frame aligned with the chordwise velocity direction
        wrench_moment = sym.Matrix([0,0,0,0,self.dM_w,0])
        velocity_frame = Transform.R_y(delta_alpha)

        self.dQ_M = (velocity_frame.ManipJacobian(p.q)).T
        self.dQ_M *= velocity_frame.InvAdjoint().T
        self.dQ_M *= wrench_moment

        self.dQ = self.dQ_L + self.dQ_M

        # generalised force
        _Q = self.dQ.integrate(int_tuple)

        tup = FwtParams.GetTuple()   
        self.dAlpha_func = sym.lambdify((tup,FwtParams.x),self.dAlpha,"numpy")

        super().__init__(p,_Q)       

    def GetAlpha(self,tup,x,t,**kwargs):
        return self.dAlpha_func(tup,x)




