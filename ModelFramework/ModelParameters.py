import sympy as sym
import sympy.physics.mechanics as me

class ModelValue:
    """
    Base class to inject a value onto sympy classes
    """
    def __init__(self,value,**kwarg):
        self.value = value
        self._dependent = False
        super().__init__(**kwarg)
    
    def __call__(self,t,x):
        return self._GetValue(t,x)

    def _GetValue(self,t,x):
        if callable(self.value):
            return self.value(t,x)
        else:
            return self.value

    def GetSub(self,t,x):
        return self(t,x)
        
class ModelSymbol(sym.Symbol,ModelValue):
    """
    Wrapper for Sympy Symbol, to inject it with a value attribute
    """
    def __init__(self,string,**kwarg):
        super().__init__(**kwarg)
    def __new__(cls,string,**kwarg):
        return super().__new__(cls,string)
    
class ModelMatrix(sym.Matrix,ModelValue):
    """
    Wrapper for Sympy Matrix, to inject it with a value attribute
    """
    def __init__(self,symbols,**kwarg):
        super().__init__(**kwarg)
    def __new__(cls,symbols,**kwargs):
        return super().__new__(cls,symbols)

class ModelExpr(sym.Symbol,ModelValue):
    def __init__(self,string,func,**kwarg):
        self.expr_func = func
        super().__init__(**kwarg)

    def _GetValue(self,t,x):
        return self.expr_func(t,x)

    def __new__(cls,string,**kwarg):
        return super().__new__(cls,string)

    def GetSub(self,t,x):
        return self.value
    

class ModelParameters:

    @classmethod
    def DynamicModel(cls,DoFs):
        model = cls()

        model.qs = DoFs
        model.q = sym.Matrix(me.dynamicsymbols(f'q:{DoFs}'))
        model.qd = sym.Matrix(me.dynamicsymbols(f'q:{DoFs}',1))
        model.qdd = sym.Matrix(me.dynamicsymbols(f'q:{DoFs}',2))

        # create state matrix
        x_ls = []
        for i in range(0,DoFs):
            x_ls.append(model.q[i])
            x_ls.append(model.qd[i])
        model.x = sym.Matrix(x_ls)
        return model      
    
    def GetTuple(self,ignore=[]):
        return tuple(var for name,var in vars(self).items() if isinstance(var,ModelValue) and name not in ignore and var not in ignore)
    
    def GetSubs(self,t,x,ignore=[]):
        sub_dependent_dict = {}
        sub_dict = {}
        # put dependent substitions in first
        for name,var in vars(self).items():
            if isinstance(var,ModelValue) and name not in ignore and var not in ignore:
                if isinstance(var,ModelMatrix):
                    for i in range(len(var)):
                        sub_dict[var[i]] = var.GetSub(t,x)[i]
                else:
                    if var._dependent:
                        sub_dependent_dict[var] = var.GetSub(t,x)
                    else:
                        sub_dict[var] = var.GetSub(t,x)
        # sub in values for dependent subsitutions
        for key,value in sub_dependent_dict.items():
            sub_dependent_dict[key] = value.subs(sub_dict)
        # return combined dictionaries
        return {**sub_dict,**sub_dependent_dict}
    
    def GetNumericTuple(self,x,t,ignore=[]):
        return tuple(var(t,x) for name,var in vars(self).items() if isinstance(var,ModelValue) and name not in ignore and var not in ignore)

    