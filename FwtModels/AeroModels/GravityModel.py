import sympy as sym
import sympy.physics.mechanics as me
import numpy as np

class GravityModel:

    def __init__(self,FwtParams,Transform,ForceVector):
        p = FwtParams
       
        #jacobian for the CoM of the FWT
        dr_idq_j= Transform.Transform_point([0,0,0]).jacobian(p.q)

        # generalised force
        self._Q = (dr_idq_j.T*ForceVector)

        #convert into a Lambda function
        self.q_func = self.GenerateLambdas(p)

    def GenerateLambdas(self,FwtParams):
        tup = FwtParams.GetTuple()
        q_func = sym.lambdify((*tup,FwtParams.x),self._Q)
        return q_func

    def __call__(self,FwtParams,x,t):
        tup = FwtParams.GetNumericTuple(x,t)
        Q_g = self.q_func(*tup,x)[:,0]
        max_len = 0
        for i in Q_g:
            if isinstance(i,float):
                len_i = 1
            else:
                len_i = len(i)  
            max_len = len_i if len_i>max_len else max_len
        if max_len == 1:
            return Q_g
        Q_g_new = np.ones((max_len,len(Q_g)))
        for i in range(0,len(Q_g)):
            Q_g_new[:,i]=Q_g[i]
        return Q_g_new.T

    def Q(self):
        return self._Q
