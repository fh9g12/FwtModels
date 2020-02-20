import sympy as sym

class FwtVariable(sym.Symbol):
    """ child class of the symbol class to pepper it with a numerical value"""
    def __init__(self,v,sStr):
        self.value = v
        super().__init__()
    def __new__(cls,v,sStr):
        return super().__new__(cls,sStr)


class FwtParameters:

    @classmethod 
    def Default2DoF(cls):
        inst = cls()
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
        return inst
    
    def GetTuple(self):
        return tuple([v for k, v in vars(self).items() if isinstance(v,sym.Symbol) ])

    def GetNumericTuple(self):
        return tuple([v.value for k, v in vars(self).items() if isinstance(v,sym.Symbol) ])
    