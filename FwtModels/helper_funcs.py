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
        tau = tau - fac*q[qi]*y_s**(i+1)
    
    z -= tau*(x-x_f)

    return z, tau

def GetVh(alpha,beta,Lambda,theta):
    ## a/c to wind transform
    Ac_V = sym.rot_axis3(beta)*sym.rot_axis2(alpha)   # transform from a/c to the velocity frame 
    ## a/c to Hinge transform
    H_Ac = sym.trigsimp(sym.rot_axis3(-Lambda)*\
                    sym.rot_axis1(theta)*  \
                    sym.rot_axis3(Lambda))        # transform from a/c to the hinge frame 


    H_V = sym.trigsimp(H_Ac*Ac_V) # transform from velocity to hinge reference frame

    # Velocity vector in velocity frame is of the form [v 0 0]
    Vv = sym.Matrix([1,0,0])
    # Transform into the hinge reference frame
    return sym.simplify(H_V * Vv)

def GetAoA(alpha,beta,Lambda,theta):
    Vh = GetVh(alpha,beta,Lambda,theta)
    return sym.atan(Vh[2]/Vh[0])
    

