import sympy as sym
from . import ExternalForce

class ContinuousGravityForce(ExternalForce):

    def __init__(self,FwtParams,Transform,ForceVector,*int_tuple):
        p = FwtParams
        # create the wrench applied at the origin of the endeffector in spetial coords
        wrench_g = sym.Matrix([*ForceVector,0,0,0])

        ## convert this to a spatial wrench

        # make a frame a end effector in body frame
        T_trans = Transform.PuesdoSpatialFrame()

        # convert wrench into this frame
        F_s = T_trans.Adjoint().T*wrench_g

        # convert into joint torques
        _dQ = sym.simplify(T_trans.ManipJacobian(p.q).T*F_s)
        _Q = _dQ.integrate(*int_tuple)

        super().__init__(p,_Q)
