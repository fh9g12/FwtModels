import sympy as sym
import sympy.physics.mechanics as me
import numpy as np

class ForcingFunction:

    def __init__(self,Forcingfunction):
        self.Forcingfunction = Forcingfunction

    def __call__(self,tup,x,t,**kwargs):
        return self.Forcingfunction(tup,x,t)

    def Q(self):
        return None
