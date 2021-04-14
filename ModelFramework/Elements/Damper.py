import sympy as sym
from .BaseElement import BaseElement

class Damper(BaseElement):
    def __init__(self,velocity,damping_constant):
        self.__c = damping_constant
        self.__z_dot = velocity
    def CalcKE(self,p):
        return 0
    def CalcPE(self,p):
        return 0
    def CalcRDF(self,p):
        return sym.Rational(1,2)*self.__c*self.__z_dot**2