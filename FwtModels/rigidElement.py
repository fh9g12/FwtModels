import sympy as sym
import numpy as np
from .base_element import BaseElement
from sympyTransforms import Vee,Wedge

def MassMatrix(m,I_xx=0,I_yy=0,I_zz=0,I_xy=0,I_xz=0,I_yz=0):
    M = sym.diag(m,m,m,I_xx,I_yy,I_zz)
    M[4,3]=M[3,4]=I_xy
    M[5,3]=M[3,5]=I_xz
    M[5,4]=M[4,5]=I_yz
    return M


class RigidElement(BaseElement):
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
        J_xyz = self.Transform.Transform_point([0,0,0]).jacobian(q)
        # get jacobian of rotations
        J_r = self.Rotations.jacobian(q)

        # construct full jacobian
        J=sym.zeros(6,len(q))
        J[:3,:] = J_xyz
        J[3:,:] = J_r

        return J
    
    
    def Jacobian2(self,q):
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