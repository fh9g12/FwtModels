import sympy as sym
import numpy as np
from .base_element import BaseElement
from sympyTransforms import Vee,Wedge

class RigidElementv2(BaseElement):
    def __init__(self,Transform,Rotations,M):
        self.Transform = Transform
        self.M_e = M
        self.Rotations = sym.Matrix([0,0,0]) if Rotations is None \
                                                else Rotations

    @classmethod
    def PointMass(cls, Transform,Rotations,m):
        return cls(Transform,Rotations,MassMatrix(m))
    
    
    def Jacobian(self,q):
        # create the jacobian for the mass
        inv = Transform.E**-1
        J = sym.Matrix([[0]*p.qs]*6)
        for i,qi in enumerate(p.q):
            J[:,i] = Vee(Transform.E.diff(qi)*inv)
        return sym.simplify(J)
    

    def CalcKE(self,p):
        # create the jacobian for the mass
        J = self.Jacobian(p.q)

        #get M in world frame
        #calculate the mass Matrix
        M = J.T*self.M_e*J

        # calculate the K.E
        T = sym.Rational(1,2)*p.qd.T*M*p.qd
        return T[0]

    def CalcPE(self,p):
        return 0