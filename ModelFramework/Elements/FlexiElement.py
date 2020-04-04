import sympy as sym
from .BaseElement import BaseElement

class FlexiElement(BaseElement):
    def __init__(self,Transform,M,x,y,z,c,s,x_f,EI,GJ,gravityPot = False):
        
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

        self._gravityPotential = gravityPot

        self.Transform = Transform
        self.dTransform = Transform.Translate(self.x,self.y,self.z)
        self.M_e = M

    def CalcKE(self, p):
        M = self.M(p)
        # calculate the K.E
        T = sym.Rational(1,2)*p.qd.T*M*p.qd
        return sym.trigsimp(T[0].integrate(self.x_integral,self.y_integral)).expand()

    
    def M(self,p):
        # create the jacobian for the mass    
        Js = self.dTransform.ManipJacobian(p.q)
        Jb = self.dTransform.InvAdjoint()*Js
        #calculate the mass Matrix in world frame
        return Jb.T*self.M_e*Jb


    def CalcElasticPE(self,p):
        # Bending Potential Energy per unit length
        v = self.dTransform.subs(self.x,self.x_f).diff(self.y,self.y).Transform_point([0,0,0])
        U_e = sym.trigsimp((v.T*v))[0]*self.EI*sym.Rational(1,2)

        # Torsional P.E per unit length
        v = self.dTransform.diff(self.x).diff(self.y).Transform_point([self.x_f,0,0])
        U_e += sym.trigsimp((v.T*v))[0]*self.GJ*sym.Rational(1,2)

        return U_e.integrate(self.y_integral)

    def CalcGravitionalPE(self, p):
        point_z = self.dTransform.Transform_point([0]*3)[2]
        dmg = point_z*self.M_e[0,0]*p.g
        return dmg.integrate(self.x_integral,self.y_integral)

    def CalcPE(self,p):
        PE = self.CalcElasticPE(p)
        PE += self.CalcGravitionalPE(p) if self._gravityPotential else sym.Integer(0)
        return PE





            
