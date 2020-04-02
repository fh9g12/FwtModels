import sympy as sym
import numpy as np
import sympy.physics.mechanics as me
import pandas as pd
from .LambdifyExtension import msub

def ShapeFunctions_BN_TM(n,m,q,y_s,x,x_f,alpha_r,factor = 1):
    # check q is the length of n+m
    if n+m != len(q):
        raise ValueError('the sum of n+m must be the same as a length of q')

    # make factor a list the size of n+m
    if isinstance(factor,int) | isinstance(factor,float):
        factor = [factor]*(n+m)
    z = sym.Integer(0)
    tau = alpha_r

    for i in range(0,n):
        z = z + q[i]*y_s**(2+i)/factor[i]
    for i in range(0,m):
        qi = i+n
        tau = tau - q[qi]*y_s**(i+1)/factor[n+i]
    
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

def LineariseMatrix(M,x,x_f):
    # reverse order of states to ensure velocities are subbed first
    x_subs = [[x[i],x_f[i],2] for i in range(len(x))]

    # get the value of M at the fixed point
    M_f = msub(M,x_subs)

    # add a gradient term for each state about the fixed point
    x_subs = {(x[i],x_f[i]) for i in range(-1,-len(x)-1,-1)}
    for i,xi in enumerate(x):
        M_f += M.diff(xi).subs(x_subs)*(xi-x_f[i])
    return M_f

def GetCruiseConditions(dof2Model,FwtParams,vs,initialGuess):
    p =FwtParams
    # data frame to append all results to
    df = pd.DataFrame(columns=['v','aoa','q'])

    # fro each velocity find the cruise state    
    for v in vs:
        p.V.value = v
        x = minimize(__ToMinimise,[initialGuess],args=(dof2Model,p)).x
        initialGuess = x
        df = df.append({'aoa' : p.alpha_r.value, 'v' : v,'q':x} , ignore_index=True)
    return df

def __ToMinimise(q,dof2Model,p):
        val = dof2Model.deriv(0,[i for i  in __CreateZeroVelStateIterable(q)],p)
        return val[1]**2 + val[3]**2

def __CreateZeroVelStateIterable(q):
    for _q in q:
        yield _q
        yield 0

def ExtractEigenValueData(evals,evecs,margin=1e-9,sortby=None):       
    # get unique eigen values
    unique = []
    unique_vecs = []
    for idx,val in enumerate(evals):
        if np.iscomplex(val):
            # check the complex conjugate is not already in the list
            if not any(np.isclose(np.conj(val),unique)):
                unique.append(val)
                unique_vecs.append(evecs[:,idx])
        else:
            # and real poles straight away
            unique.append(val)
            unique_vecs.append(evecs[:,idx].tolist())
    df_v = pd.DataFrame()
    
    df_v['cn'] = unique
    df_v['Real'] = np.real(unique)
    df_v['Imag'] = np.imag(unique)
    df_v['Frequency'] = np.where(np.iscomplex(unique),np.abs(unique)/(2*np.pi),0)
    df_v['Damping'] = np.where(np.iscomplex(unique),np.cos(np.angle(unique)),np.NaN)
    df_v['Stable'] = np.max(df_v['Real'])<=margin
    df_v['Eigen Vector'] = unique_vecs

    if sortby is not None:
        df_v = df_v.sort_values(by=sortby)  
    df_v['Mode'] = range(0,len(unique))
    return df_v
