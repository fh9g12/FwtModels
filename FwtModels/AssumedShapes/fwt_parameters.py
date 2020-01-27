import sympy as sym

class FwtVariable(sym.Symbol):
    """ child class of the symbol class to pepper it with a numerical value"""
    def __init__(self,v,sStr):
        self.value = v
        super().__init__()
    def __new__(cls,v,sStr):
        return super().__new__(cls,sStr)

class FwtParameters:
    m_w: FwtVariable = FwtVariable(0,'m_w')
    m_t: FwtVariable = FwtVariable(0,'m_t')
    x_f: FwtVariable = FwtVariable(0,'x_f')
    s_w: FwtVariable = FwtVariable(0,'s_w')
    s_t: FwtVariable = FwtVariable(0,'s_t')
    c: FwtVariable = FwtVariable(0,'c')
    Lambda: FwtVariable = FwtVariable(0,'Lambda')
    EI: FwtVariable = FwtVariable(0,'EI')
    GJ: FwtVariable = FwtVariable(0,'GJ')
    k_theta : FwtVariable = FwtVariable(0,'k_theta')
    rho: FwtVariable = FwtVariable(0,'rho')
    V: FwtVariable = FwtVariable(0,'V')
    a_w : FwtVariable = FwtVariable(0,'a_w')
    a_t : FwtVariable = FwtVariable(0,'a_t')
    alpha_0 : FwtVariable = FwtVariable(0,'alpha_0')
    e : FwtVariable = FwtVariable(0,'e')
    Malphadot : FwtVariable = FwtVariable(0,'M_alphadot')
    g : FwtVariable = FwtVariable(0,'g')
        
    
    def __init__(self,m_w = 0,m_t = 0,x_f=0,s_w=0,s_t=0,c=0,Lambda=0,EI=0,GJ=0,
        k_theta=0,rho=0,V=0,a_w=0,a_t=0,alpha_0=0,e=0,Malphadot=0,g = 9.81):
        self.m_w.value = m_w
        self.m_t.value = m_t
        self.x_f.value = x_f
        self.s_w.value = s_w
        self.s_t.value = s_t
        self.c.value = c
        self.Lambda.value = Lambda
        self.EI.value = EI
        self.GJ.value = GJ
        self.k_theta.value = k_theta
        self.rho.value = rho
        self.V.value = V
        self.a_w.value = a_w
        self.a_t.value = a_t
        self.alpha_0.value = alpha_0
        self.e.value = e
        self.Malphadot.value = Malphadot
        self.g.value = g
    
    def GetTuple(self):
        return (self.m_w,self.m_t,self.x_f,self.s_w,self.s_t,self.c,
        self.Lambda,self.EI,self.GJ,self.k_theta,self.rho,self.V,
        self.a_w,self.a_t,self.alpha_0,self.e,self.Malphadot,self.g)
    