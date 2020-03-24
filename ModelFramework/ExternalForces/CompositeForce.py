import sympy as sym
import numpy as np
from . import ExternalForce

class CompositeForce(ExternalForce):

    def __init__(self,p,forces):
        self.forces = forces
        self.__qs = p.qs

    def __call__(self,tup,x,t):
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

    def subs(self,p,*args):
        new_forces = []
        for force in self.forces:
            new_forces.append(force.subs(p,*args))
        return CompositeForce(p,new_forces)

    def gensource(self,name = 'Q'):
        # add each force to the string
        full_str = ''
        for force,i in self.forces:
            full_str += force.gensource(f'Q_{i}')
        # Add the main force
        lines = []
        lines.append(f'def {name}(tup,x,t):')
        lines.append(f'\tval = zeros(({self.__qs},1))')
        for i in len(self.forces):
            lines.append(f'\tval += Q_{i}(tup,x,t)')
        lines.append('\treturn val')
        return '\n'.join(lines)+'\n'
