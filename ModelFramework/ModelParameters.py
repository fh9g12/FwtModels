import sympy as sym
import sympy.physics.mechanics as me

class ModelVariable(sym.Symbol):
    """ child class of the symbol class to pepper it with a numerical value"""
    def __init__(self,v,sStr):
        self.value = v
        super().__init__()
    def __new__(cls,v,sStr):
        return super().__new__(cls,sStr)


class ModelParameters:

    @classmethod
    def DynamicModel(cls,DoFs):
        model = cls()

        model.qs = qs
        model.q = sym.Matrix(me.dynamicsymbols(f'q:{qs}'))
        model.qd = sym.Matrix(me.dynamicsymbols(f'q:{qs}',1))
        model.qdd = sym.Matrix(me.dynamicsymbols(f'q:{qs}',2))

        # create state matrix
        x_ls = []
        for i in range(0,qs):
            x_ls.append(model.q[i])
            x_ls.append(model.qd[i])
        model.x = sym.Matrix(x_ls)        
    
    def GetTuple(self,ignore=[]):
        return tuple([v for k,v in vars(self).items() if isinstance(v,ModelVariable) and k not in ignore and v not in ignore ])
    
    def GetSubs(self,ignore=[]):
        return {v:v.value for k,v in vars(self).items() if isinstance(v,ModelVariable) and k not in ignore and v not in ignore}

    def GetNumericTuple(self,x,t,ignore=[]):
        vals = []
        for k, v in vars(self).items():
            if isinstance(v,ModelVariable) and (k not in ignore or v not in ignore):
                vals.append(v.value(t,x) if callable(v.value) else v.value)
        return tuple(vals)
        


    