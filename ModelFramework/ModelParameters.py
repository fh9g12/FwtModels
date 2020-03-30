import sympy as sym
import sympy.physics.mechanics as me

class ModelValue:
    """
    Base class to inject a value onto sympy classes
    """
    def __init__(self,value,**kwarg):
        self.value = value
        super().__init__(**kwarg)
    
    #def __call__(self,t,x):
    #    return self.value

    #def GetValue
        
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

#class ModelExpr(ModelValue):
#    def __init__(self,func,**kwarg):
#        self.expr_func = func
#    def __call__(self,t,x):
#        return self.expr_func(t,x)
    

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
        return tuple([v for k,v in vars(self).items() if isinstance(v,ModelValue) and k not in ignore and v not in ignore ])
    
    def GetSubs(self,ignore=[]):
        return {v:v.value for k,v in vars(self).items() if isinstance(v,ModelVaModelValueriable) and k not in ignore and v not in ignore}

    def GetNumericTuple(self,x,t,ignore=[]):
        vals = []
        for k, v in vars(self).items():
            if isinstance(v,ModelValue) and (k not in ignore or v not in ignore):
                vals.append(v.value(t,x) if callable(v.value) else v.value)
        return tuple(vals)
        


    