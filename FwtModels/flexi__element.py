import sympy as sym
import numpy as np
from .base_element import BaseElement

class FlexiElement(BaseElement):
    def __init__(self,Transform,Rotations,M,x,y,z,c,s,x_f,EI,GJ):
        
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
        self.Rotations = sym.Matrix([0,0,0]) if Rotations is None \
                                                else Rotations

    def Jacobian(self,q):
        # create the jacobian for the mass
        J_xyz = self.Transform.Transform_point([0,0,0]).jacobian(q)
        # get jacobian of rotations
        J_r = self.Rotations.jacobian(q)
        # construct full jacobian
        J=sym.zeros(6,len(q))
        J[:3,:] = J_xyz
        J[3:,:] = J_r

        return J

    def CalcKE(self, p):
        # create the jacobian for the mass
        J = self.Jacobian(p.q)

        #get M in world frame
        #calculate the mass Matrix
        M = J.T*self.M_e*J

        # calculate the K.E
        T = sym.Rational(1,2)*p.qd.T*M*p.qd
        return sym.simplify(T[0].integrate((self.x,-self.c,0),(self.y,0,self.s)))

    def CalcPE(self,p):
        #calc pot energy per unit length
        U_e = sym.Rational(1,2)*(self.z.subs(self.x,self.x_f).diff(self.y,self.y)**2*self.EI)
        U_e = U_e + sym.Rational(1,2)*(self.Rotations[1].diff(self.y)**2*self.GJ)

        return U_e.integrate((self.y,0,self.s))


            
