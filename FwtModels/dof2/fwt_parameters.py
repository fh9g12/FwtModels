import sympy as sym

class FwtVariable(sym.Symbol):
    """ child class of the symbol class to pepper it with a numerical value"""
    def __init__(self,v,sStr):
        self.value = v
        super().__init__()
    def __new__(cls,v,sStr):
        return super().__new__(cls,sStr)


class FwtParameters:
    m: FwtVariable = FwtVariable(0,'m') # mass of FWT
    l: FwtVariable = FwtVariable(0,'l') # dist from hinge to CoM
    s: FwtVariable = FwtVariable(0,'s') # span
    c: FwtVariable = FwtVariable(0,'c') # chord
    k: FwtVariable = FwtVariable(0,'k') # spring constant
    g : FwtVariable = FwtVariable(0,'g') # gravity
    Lambda: FwtVariable = FwtVariable(0,'Lambda') # flare angle
    rho: FwtVariable = FwtVariable(0,'rho') # density
    V: FwtVariable = FwtVariable(0,'V') # velocity
    a_t : FwtVariable = FwtVariable(0,'a_t') # C_L of FWT
      
    def __init__(self):
        self.params = []
        self.params.append(self.m)
        self.params.append(self.l)
        self.params.append(self.s)
        self.params.append(self.c)
        self.params.append(self.k)
        self.params.append(self.g)
        self.params.append(self.Lambda)
        self.params.append(self.rho)
        self.params.append(self.V)
        self.params.append(self.a_t)

    def AddParam(self,value,symString):
        param = FwtVariable(value,symString)
        self.params.append(param)
        return param
    
    def GetTuple(self):
        return tuple(s for s in self.params)

    def GetNumericTuple(self):
        return tuple(s.value for s in self.params)
    