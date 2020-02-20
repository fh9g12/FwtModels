import sympy as sym
import sympy.physics.mechanics as me

class SimpleAeroModel:

    def __init__(self,FwtParams,Transform,q,qd,at_mode=0):
        p = FwtParams
        ## force per unit length will follow theredosons unsteady theory

        # symbol for y of FWT
        y_f = sym.Symbol('y_f')

        # aero focre equation
        _alpha_r, V_t = me.dynamicsymbols('alpha_r V_t')
        half = sym.Rational(1,2)

        self.dAlpha = sym.atan(sym.sin(q[0]-sym.pi*half)*sym.sin(p.Lambda))
        self.dAlpha = self.dAlpha - ((p.s*half*qd[0]+qd[0]*y_f))/V_t
        self.dAlpha = self.dAlpha - (qd[1]*sym.sin(q[0]))/V_t
        self.dAlpha = self.dAlpha + _alpha_r*sym.cos(q[0]-sym.pi*half)
        self.dL_w = half*p.rho*V_t**2*p.c*p.a_t*self.dAlpha
        self.dL_w = self.dL_w*sym.cos(self.dAlpha.subs(y_f,0))
        #jacobian per unit length of FWT
        dr_idq_j= Transform.Transform_point([0,y_f,0]).jacobian(q)

        # force Rotated into task space
        self.dL_wi = Transform.T[:3,:3]*sym.Matrix([0,0,self.dL_w])

        # generalised force per unit length
        self.dQ = dr_idq_j.T*self.dL_wi

        # generalised force
        self.Q = self.dQ.integrate((y_f,p.s*-half,p.s*half))

        #convert into a Lambda function
        tup = FwtParams.GetTuple()
        self.func = sym.lambdify((*tup,_alpha_r, V_t,q,qd),self.Q)

        self._v_func = lambda t: 20
        self._ar_func = lambda t: 0

    def SetFuncs(self,FwtParams,V_Fuc,a_r_func):
        self._v_func = V_Fuc
        self._ar_func = a_r_func


    def Calc(self,FwtParams,q,qd,t,v_0=20):
        tup = FwtParams.GetNumericTuple()
        return self.func(*tup,self._ar_func(t),self._v_func(t),q,qd)



