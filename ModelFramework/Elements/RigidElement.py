import sympy as sym
from .BaseElement import BaseElement
from .MassMatrix import MassMatrix

class RigidElement(BaseElement):
    def __init__(self,Transform,M,gravityPotential=False):
        self._gravityPotential = gravityPotential
        self.Transform = Transform
        self.M_e = M

    @classmethod
    def PointMass(cls, Transform,m,gravityPotential=False):
        return cls(Transform,MassMatrix(m),gravityPotential)    
    
    def CalcKE(self,p):
        M = self.M(p)   

        # calculate the K.E
        T = sym.Rational(1,2)*p.qd.T*M*p.qd
        return sym.trigsimp(T[0])

    def M(self,p):
        # create the jacobian for the mass
        Js = self.Transform.ManipJacobian(p.q)
        Jb = self.Transform.InvAdjoint()*Js
        #get M in world frame
        #calculate the mass Matrix
        return Jb.T*self.M_e*Jb


    def CalcPE(self,p):
        if self._gravityPotential:
            #return 0
            p = self.Transform.Transform_point([0,0,0])
            return p[2]*self.M_e[0,0]*9.81
        else:
            return 0