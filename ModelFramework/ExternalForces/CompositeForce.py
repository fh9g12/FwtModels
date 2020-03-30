import sympy as sym
import numpy as np
from . import ExternalForce

class CompositeForce(ExternalForce):

    def __init__(self,p,forces):
        self.forces = forces
        self.__qs = p.qs

    def __call__(self,tup,x,t):
        return sum((force(tup,x,t) for force in self.forces))

    def Q(self):
        val = sym.Matrix([0]*self.__qs)
        for force in self.forces:
            new_val = force.Q()
            if new_val is not None:        
                val += new_val
        return val

    def subs(self,p,*args):
        return CompositeForce(p,[force.subs(p,*args) for force in self.forces])

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

    def linearise(self,p):
        return CompositeForce(p,[force.linearise(p) for force in self.forces])
