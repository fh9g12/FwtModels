import sympy as sym
import numpy as np
import sympy.physics.mechanics as me

class CompositeForce:

    def __init__(self,p,forces):
        self.forces = forces
        self.__qs = p.qs

    def __call__(self,tup,x,t,**kwargs):
        result = np.array([[0] for i in range(0,self.__qs)])
        for i in range(0,len(self.forces)):   
            result = result + self.forces[i](tup,x,t) 
        return result

    def Q(self):
        val = sym.Matrix([0]*self.__qs)
        for i in range(0,len(self.forces)):
            new_val = self.forces[i].Q()
            if new_val is not None:        
                val += new_val
        return val