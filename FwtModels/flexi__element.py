import sympy as sym
import numpy as np
from .base_element import BaseElement
from sympyTransforms import Vee,Wedge

class FlexiElement(BaseElement):
    def __init__(self,Transform,M,x,y,z,c,s,x_f,EI,GJ):
        
        self.x = x
        self.y = y
        self.z = z
        self.c = c
        self.s = s
        self.EI = EI
        self.GJ = GJ
        self.x_f = x_f

        self.Transform = Transform
        self.M_e = M


    #def Jacobian(self,q):
    #    # create the jacobian for the mass
    #    inv = self.Transform.E**-1
    #    J = sym.zeros(6,len(q))
    #    for i,qi in enumerate(q):
    #        J[:,i] = Vee(self.Transform.E.diff(qi)*inv)
    #    return sym.simplify(J)

    def CalcKE(self, p):
        # create the jacobian for the mass
        
        J = self.Transform.Translate(self.x,self.y,self.z).ManipJacobian(p.q)

        #get M in world frame
        #calculate the mass Matrix
        M = J.T*self.M_e*J

        # calculate the K.E
        T = sym.Rational(1,2)*p.qd.T*M*p.qd
        return sym.simplify(T[0].integrate((self.x,0,self.c),(self.y,0,self.s)))

    def CalcPE(self,p):
        #first derivative
        Trans = self.Transform.Translate(self.x_f,self.y,self.z.subs(self.x,self.x+self.x_f))

        # Bending Potential Energy per unit length
        v = Trans.diff(self.y).diff(self.y).Transform_point([0,0,0])
        U_e = sym.trigsimp((v.T*v))[0]*self.EI*sym.Rational(1,2)

        # Torsional P.E per unit length
        v = Trans.diff(self.x).diff(self.y).Transform_point([0,0,0])
        U_e += sym.trigsimp((v.T*v))[0]*self.GJ*sym.Rational(1,2)

        return U_e.integrate((self.y,0,self.s))


            
