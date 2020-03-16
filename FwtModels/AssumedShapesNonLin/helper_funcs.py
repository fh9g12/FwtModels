import sympy as sym
import numpy as np
import sympy.physics.mechanics as me

def ShapeFunctions_BN_TM(n,m,q,y_s,x,x_f,alpha_r,factor = 1000):
    # check q is the length of n+m
    if n+m != len(q):
        raise ValueError('the sum of n+m must be the same as a length of q')
    fac = sym.Rational(1,factor)

    z = sym.Integer(0)
    tau = alpha_r

    for i in range(0,n):
        z = z + fac*q[i]*y_s**(2+i)
    for i in range(0,m):
        qi = i+n
        tau = tau + fac*q[qi]*y_s**(i+1)
    
    z = z - tau*(x-x_f)

    return z, tau
    

