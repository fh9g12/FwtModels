import sympy as sym
from inspect import getsource
from ..LambdifyExtension import msub
from ..helper_funcs import LineariseMatrix


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
        return ExternalForce(p,LineariseMatrix(self.Q(),p.x,p.fp))





