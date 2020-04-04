import sympy as sym
from .BaseElement import BaseElement

class Spring(BaseElement):
    def __init__(self,deflection,spring_constant):
        self.__k = spring_constant
        self.__z = deflection
    def CalcKE(self,p):
        return 0
    def CalcPE(self,p):
        return sym.Rational(1,2)*self.__k*self.__z**2
        
