import sympy as sym
from .BaseElement import BaseElement

class FlexiElement(BaseElement):
    def __init__(self,Transform,M,x,y,z,c,s,x_f,EI,GJ,c_negative=False,s_negative=False):
        
        self.x = x
        self.y = y
        self.z = z
        self.c = c
        self.s = s
        self.EI = EI
        self.GJ = GJ
        self.x_f = x_f

        self.y_integral = (self.y,0,self.s)
        self.x_integral = (self.x,0,self.c) 

        self.Transform = Transform
        self.M_e = M

    def CalcKE(self, p):
        M = self.M(p)
        # calculate the K.E
        T = sym.Rational(1,2)*p.qd.T*M*p.qd
        return sym.simplify(T[0].integrate(self.x_integral,self.y_integral))

    def M(self,p):
        # create the jacobian for the mass    
        T = self.Transform.Translate(self.x,self.y,self.z)
        Js = T.ManipJacobian(p.q)
        Jb = T.InvAdjoint()*Js
        #get M in world frame
        #calculate the mass Matrix
        return Jb.T*self.M_e*Jb


    def CalcPE(self,p):
        #first derivative
        Trans = self.Transform.Translate(self.x_f,self.y,self.z.subs(self.x,self.x+self.x_f))

        # Bending Potential Energy per unit length
        v = Trans.subs(self.x,0).diff(self.y,self.y).Transform_point([0,0,0])
        U_e = sym.trigsimp((v.T*v))[0]*self.EI*sym.Rational(1,2)

        # Torsional P.E per unit length
        v = Trans.diff(self.x).diff(self.y).Transform_point([0,0,0])
        U_e += sym.trigsimp((v.T*v))[0]*self.GJ*sym.Rational(1,2)

        return U_e.integrate(self.y_integral)


            
