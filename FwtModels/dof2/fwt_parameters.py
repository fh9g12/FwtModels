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

    @classmethod 
    def Default2DoF(cls):
        inst = cls(2) # 2DoF System
        inst.m: FwtVariable = FwtVariable(0,'m') # mass of FWT
        inst.l: FwtVariable = FwtVariable(0,'l') # dist from hinge to CoM
        inst.s: FwtVariable = FwtVariable(0,'s') # span
        inst.c: FwtVariable = FwtVariable(0,'c') # chord
        inst.k: FwtVariable = FwtVariable(0,'k') # spring constant
        inst.g : FwtVariable = FwtVariable(0,'g') # gravity
        inst.Lambda: FwtVariable = FwtVariable(0,'Lambda') # flare angle
        inst.rho: FwtVariable = FwtVariable(0,'rho') # density
        inst.V: FwtVariable = FwtVariable(0,'V') # velocity
        inst.a_t : FwtVariable = FwtVariable(0,'a_t') # C_L of FWT
        inst.alpha_r : FwtVariable = FwtVariable(0,'alpha_r') # C_L of FWT
        inst.q : sym.Matrix(me.dynamicsymbols(f'q:{2}'))
        return inst
    
    def GetTuple(self):
        return tuple([v for k, v in vars(self).items() if isinstance(v,sym.Symbol) ])

    def GetNumericTuple(self,x,t):
        vals = []
        for _, v in vars(self).items():
            if isinstance(v,sym.Symbol):
                vals.append(v.value(t,x) if callable(v.value) else v.value)
        return tuple(vals)
    