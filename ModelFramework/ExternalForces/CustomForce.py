from . import ExternalForce

class CustomForce(ExternalForce):

    def __init__(self,Params,Forcingfunction):
        self.q_func = Forcingfunction
        self.__qs = Params.qs

    def Q(self):
        return None
    
    def subs(self,p,*args):
        return self

    def gensource(self,name = 'Q'):
        # Add the main force
        lines = []
        lines.append(f'def {name}(tup,x,t):')
        lines.append(f'\treturn zeros(({self.__qs},1))')
        return '\n'.join(lines)+'\n'

    def linearise(self,p):
        return self

        
