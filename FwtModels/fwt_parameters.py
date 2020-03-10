import sympy as sym
import sympy.physics.mechanics as me

class FwtVariable(sym.Symbol):
    """ child class of the symbol class to pepper it with a numerical value"""
    def __init__(self,v,sStr):
        self.value = v
        super().__init__()
    def __new__(cls,v,sStr):
        return super().__new__(cls,sStr)


class FwtParameters:
    
    def __init__(self,qs):
        self.qs = qs
        self.q = sym.Matrix(me.dynamicsymbols(f'q:{qs}'))
        self.qd = sym.Matrix(me.dynamicsymbols(f'q:{qs}',1))
        self.qdd = sym.Matrix(me.dynamicsymbols(f'q:{qs}',2))

        # create state matrix
        x_ls = []
        for i in range(0,qs):
            x_ls.append(self.q[i])
            x_ls.append(self.qd[i])
        self.x = sym.Matrix(x_ls)
    
    def GetTuple(self,ignore=[]):
        return tuple([v for k,v in vars(self).items() if isinstance(v,sym.Symbol) and k not in ignore and v not in ignore ])

    def GetNumericTuple(self,x,t,ignore=[]):
        vals = []
        for k, v in vars(self).items():
            if isinstance(v,sym.Symbol) and (k not in ignore or v not in ignore):
                vals.append(v.value(t,x) if callable(v.value) else v.value)
        return tuple(vals)
        


    