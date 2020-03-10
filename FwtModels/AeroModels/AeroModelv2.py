import sympy as sym
import sympy.physics.mechanics as me
from scipy import integrate

def FwtAoA(FwtParams,foldAngle):
    p = FwtParams
    # get velocity vector in FWT frame
    v_x = p.V * (sym.sin(p.Lambda)**2*sym.cos(p.alpha_r)*sym.cos(foldAngle) - sym.sin(p.Lambda)**2*sym.cos(p.alpha_r) - sym.sin(p.Lambda)*sym.sin(p.alpha_r)*sym.sin(foldAngle) + sym.cos (p.alpha_r))
    v_z = p.V * (sym.sin(p.Lambda)*sym.sin(foldAngle)*sym.cos(p.alpha_r) + sym.sin(p.alpha_r)*sym.cos(foldAngle))
    return sym.atan(v_z/v_x)


class AeroModelv2:

    def __init__(self,FwtParams,Transform,C_L,int_tuple,AoA_expr=None):
        p = FwtParams
        ## force per unit length will following theredosons pseado-steady theory
        AoA_expr = AoA_expr if AoA_expr is not None else p.alpha_r

        # add z velocity due to motion
        v_z_eff = sym.simplify(Transform.BodyVelocity()[2])

        # combine to get effective AoA
        self.dAlpha = AoA_expr + v_z_eff/p.V

        # Calculate the lift force
        self.dL_w = -sym.Rational(1,2)*p.rho*p.V**2*p.c*C_L*self.dAlpha

        # Rotate to be perp with Hinge
        #self.dL_w = self.dL_w*sym.cos(self.dAlpha)

        #jacobian per unit length of FWT
        dr_idq_j= Transform.Transform_point([0,0,0]).jacobian(p.q)

        # force Rotated into task space
        self.dL_wi = Transform.R*sym.Matrix([0,0,self.dL_w])

        # generalised force per unit length
        self.dQ = (dr_idq_j.T*self.dL_wi)

        # generalised force
        self._Q = self.dQ.integrate(int_tuple)

        #convert into a Lambda function
        self.q_func,self.dAlpha_func = self.GenerateLambdas(p)

    def GenerateLambdas(self,FwtParams):
        tup = FwtParams.GetTuple()
        q_func = sym.lambdify((*tup,FwtParams.x),self._Q)
        dAlpha_func = sym.lambdify((*tup,FwtParams.x),self.dAlpha)
        return q_func,dAlpha_func

    def GetAlpha(self,tup,x,t,**kwargs):
        return self.dAlpha_func(*tup,x)

    def __call__(self,tup,x,t,**kwargs):
        return self.q_func(*tup,x)[:,0]

    def Q(self):
        return self._Q




