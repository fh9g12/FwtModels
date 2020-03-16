import sympy as sym
import numpy as np
from .base_element import BaseElement
from sympyTransforms import Vee,Wedge
from .mass_matrix import MassMatrix

class RigidElement(BaseElement):
    def __init__(self,Transform,M):
        self.Transform = Transform
        self.M_e = M

    @classmethod
    def PointMass(cls, Transform,m):
        return cls(Transform,MassMatrix(m))    
    
    def CalcKE(self,p):
        # create the jacobian for the mass
        J = self.Transform.ManipJacobian(p.q)

        #get M in world frame
        #calculate the mass Matrix
        M = J.T*self.M_e*J

        # calculate the K.E
        T = sym.Rational(1,2)*p.qd.T*M*p.qd
        return T[0]

    def CalcPE(self,p):
        return 0