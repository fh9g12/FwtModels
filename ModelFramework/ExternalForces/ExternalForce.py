import sympy as sym

class ExternalForce:

    def __init__(self,FwtParams,Q):
        self._Q = Q
        tup = FwtParams.GetTuple()     
        self.q_func = sym.lambdify((tup,FwtParams.x),self._Q,"numpy")

    def __call__(self,tup,x,t):
        return self.q_func(tup,x)

    def Q(self):
        return self._Q

    def subs(self,p,*args):
        return ExternalForce(p,self._Q.subs(*args))

