import sympy as sym
from . import ExternalForce

class GravityForce(ExternalForce):

    def __init__(self,FwtParams,Transform,ForceVector):
        p = FwtParams
        # create the wrench applied at the origin of the endeffector in spetial coords
        wrench_g = sym.Matrix([ForceVector[0],ForceVector[1],ForceVector[2],0,0,0])

        ## convert this to a spatial wrench

        # make a frame a end effector in body frame
        tup = tuple(Transform.t)
        T_trans = Transform.PuesdoSpatialFrame()

        # convert wrench into this frame
        F_s = T_trans.Adjoint().T*wrench_g

        # convert into joint torques
        _Q = sym.simplify(T_trans.ManipJacobian(p.q).T*F_s)

        super().__init__(p,_Q)
