import sympy as sym
from inspect import getsource
from ..LambdifyExtension import msub
from ..helper_funcs import LineariseMatrix


class ExternalForce:

    def __init__(self,Q = None):
        self._Q = Q

    def __mul__(self,other):
        return ExternalForce(self._Q*other)

    def Q(self):
        return self._Q

    def subs(self,*args):
        return ExternalForce(self._Q.subs(*args))

    def integrate(self,*args):
        return ExternalForce(self._Q.integrate(*args))

    def linearise(self,x,x_f):
        return ExternalForce(LineariseMatrix(self.Q(),x,x_f))

    def lambdify(self,params):
        if self._Q is None:
            return None
        return sym.lambdify(params,self._Q ,"numpy")




    





