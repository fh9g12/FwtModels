import sympy as sym
import sympy.physics.mechanics as me

class SimpleAeroModel:

    def __init__(self,FwtParams,Transform,at_mode=0):
        p = FwtParams
        ## force per unit length will follow theredosons unsteady theory

        # symbol for y of FWT
        self.y_t = sym.Symbol('y_t')

        # set a_t mode
        if at_mode == 1:
            C_L_expr = p.a_t - p.a_t*self.y_t/p.s
        else:
            C_L_expr = p.a_t

        # aero focre equation
        self._alpha_r, self.V_t = me.dynamicsymbols('alpha_r V_t')
        half = sym.Rational(1,2)

        self.dAlpha = sym.atan(sym.sin(p.q[0]-sym.pi*half)*sym.sin(p.Lambda))

        thet = p.q[0]-sym.pi*sym.Rational(1,2)
        stheta = sym.sin(thet)
        ctheta = sym.cos(thet)
        sLambda = sym.sin(p.Lambda)
        n = self._alpha_r*ctheta+sLambda*stheta
        d = self._alpha_r*sLambda*stheta-sLambda**2*ctheta+sLambda**2-1
        self.dAlpha = sym.atan(-n/d)

        self.dAlpha = self.dAlpha - p.qd[0]*self.y_t/self.V_t
        self.dAlpha = self.dAlpha - (p.qd[1]*sym.sin(p.q[0]))/self.V_t
        self.dL_w = half*p.rho*self.V_t**2*p.c*C_L_expr*self.dAlpha

        #self.dL_w = self.dL_w*sym.cos(self.dAlpha.subs(y_f,0))
        #jacobian per unit length of FWT
        dr_idq_j= Transform.Transform_point([0,self.y_t,0]).jacobian(p.q)

        # force Rotated into task space
        self.dL_wi = Transform.T[:3,:3]*sym.Matrix([0,0,self.dL_w])

        # generalised force per unit length
        self.dQ = (dr_idq_j.T*self.dL_wi)

        # generalised force
        self.Q = self.dQ.integrate((self.y_t,0,p.s))

        #convert into a Lambda function
        self.q_func,self.dAlpha_func = self.GenerateLambdas(p)
        
        self._v_func = lambda t: 20
        self._ar_func = lambda t: 0

    def GenerateLambdas(self,FwtParams):
        tup = FwtParams.GetTuple()
        q_func = sym.lambdify((*tup,self._alpha_r, self.V_t,FwtParams.x),self.Q)
        dAlpha_func = sym.lambdify((*tup,self._alpha_r, self.V_t,self.y_t,FwtParams.x),self.dAlpha)
        return q_func,dAlpha_func

    def SetFuncs(self,FwtParams,V_Fuc,a_r_func):
        self._v_func = V_Fuc
        self._ar_func = a_r_func

    def GetAlpha(self,FwtParams,x,y_t,t):
        tup = FwtParams.GetNumericTuple()
        return self.dAlpha_func(*tup,self._ar_func(t),self._v_func(t),y_t,x)

    def __call__(self,FwtParams,x,t):
        tup = FwtParams.GetNumericTuple()
        return self.q_func(*tup,self._ar_func(t),self._v_func(t),x)[:,0]




