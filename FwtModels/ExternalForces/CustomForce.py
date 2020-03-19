import sympy as sym
import sympy.physics.mechanics as me
import numpy as np
from . import ExternalForce

class ForcingFunction(ExternalForce):

    def __init__(self,Forcingfunction):
        self.Forcingfunction = Forcingfunction

    def __call__(self,tup,x,t):
        return self.Forcingfunction(tup,x,t)

    def Q(self):
        return None
    
    def subs(self,p,*args):
        return self
