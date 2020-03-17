import sympy as sym
import sympy.physics.mechanics as me

class AeroModelv1:

    def __init__(self,FwtParams,Transform,alpha,C_L_expr,int_tuple):
        p = FwtParams
        ## force per unit length will follow theredosons unsteady theory
        self.dAlpha = alpha

        v_z_eff = sym.simplify(Transform.BodyVelocity()[2])

        self.dAlpha += v_z_eff/p.V

        # Calculate the lift force
        self.dL_w = -sym.Rational(1,2)*p.rho*p.V**2*p.c*C_L_expr*self.dAlpha

        # Calulate the pitching Moment
        self.dM_w = sym.Rational(1,2)*p.rho*p.V**2*p.c**2
        self.dM_w *= (p.e*C_L_expr*self.dAlpha )

        self.dQ = (Transform.ManipJacobian(p.q)).T
        self.dQ *= Transform.InvAdjoint().T
        self.dQ = sym.trigsimp(self.dQ) # generally some sin^2 + cos^2 in here
        self.dQ *= sym.Matrix([0,0,self.dL_w,0,self.dM_w,0])

        # generalised force
        self._Q = self.dQ.integrate(int_tuple)

        #convert into a Lambda function
        self.q_func,self.dAlpha_func = self.GenerateLambdas(p)

    def GenerateLambdas(self,FwtParams):
        tup = FwtParams.GetTuple()
        q_func = sym.lambdify((*tup,FwtParams.x),self._Q)
        dAlpha_func = sym.lambdify((*tup,FwtParams.x),self.dAlpha)
        return q_func,dAlpha_func

    def GetAlpha(self,FwtParams,x,y_t,t):
        tup = FwtParams.GetNumericTuple(x,t)
        return self.dAlpha_func(*tup,y_t,x)

    def __call__(self,tup,x,t,**kwrags):
        #tup = FwtParams.GetNumericTuple(x,t)
        return self.q_func(*tup,x)

    def Q(self):
        return self._Q




