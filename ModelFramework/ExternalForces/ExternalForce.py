import sympy as sym
from inspect import getsource
from ..LambdifyExtension import msub


class ExternalForce:

    def __init__(self,FwtParams,Q):
        from sympy.abc import t     
        self._Q = Q
        tup = FwtParams.GetTuple()
        self.q_func = sym.lambdify((tup,FwtParams.x,t),self._Q,"numpy")

    def __call__(self,tup,x,t):
        return self.q_func(tup,x,t)

    def Q(self):
        return self._Q

    def subs(self,p,*args):
        return ExternalForce(p,self._Q.subs(*args))

    def gensource(self,name = 'Q'):
        source_str = getsource(self.q_func)
        return source_str.replace('_lambdifygenerated',name)

    def linearise(self,p):
        x_subs = {(p.x[i],p.fp[i]) for i in range(-1,-len(p.x)-1,-1)}
        Q = self.Q()
        Q_p = Q.subs(x_subs)
        for i,x in enumerate(p.x):
            Q_p += Q.diff(x).subs(x_subs)*(x-p.fp[i])
        return ExternalForce(p,Q_p)





