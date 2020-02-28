import sympy as sym
import sympy.physics.mechanics as me

class CompositeForce:

    def __init__(self,forces):
        self.forces = forces

    def __call__(self,FwtParams,x,t):
        result = self.forces[0](FwtParams,x,t)
        for i in range(1,len(self.forces)):   
            result = result +self.forces[i](FwtParams,x,t) 
        return result

    def Q(self):
        val = self.forces[0].Q()
        for i in range(1,len(self.forces)):           
            val = val + self.forces[i].Q()
        return val