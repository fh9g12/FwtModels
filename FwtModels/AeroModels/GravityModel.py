import sympy as sym
import sympy.physics.mechanics as me
import numpy as np

import sys, os
sys.path.insert(1, os.path.join(sys.path[0], '../..'))
import sympyTransforms as symt

class GravityModel:

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
        self._Q = sym.simplify(T_trans.ManipJacobian(p.q).T*F_s)
     
        self.q_func = self.GenerateLambdas(p)

    def GenerateLambdas(self,FwtParams):
        tup = FwtParams.GetTuple()
        q_func = sym.lambdify((tup,FwtParams.x),self._Q,"numpy")
        return q_func

    def __call__(self,tup,x,t,**kwargs):
        return self.q_func(tup,x)

    def Q(self):
        return self._Q
