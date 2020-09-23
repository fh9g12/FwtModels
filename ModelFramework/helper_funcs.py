import sympy as sym
import numpy as np
import sympy.physics.mechanics as me
from sympy.abc import t
import pandas as pd

def LineariseMatrix(M,x,x_f):
    # reverse order of states to ensure velocities are subbed first
    x_subs = {x[i]:x_f[i] for i in range(len(x))}

    # get the value of M at the fixed point
    M_f = me.msubs(M,x_subs)

    # add a gradient term for each state about the fixed point
    
    for i,xi in enumerate(x):
        M_f += me.msubs(M.diff(xi),x_subs)*(xi-x_f[i])
    return M_f

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

def ExtractEigenValueData_list(evals,evecs,margin=1e-9,sortby=None):       
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

    # Generate data
    real = np.real(unique)
    imag = np.imag(unique)
    F = np.where(np.iscomplex(unique),np.abs(unique)/(2*np.pi),0)
    D = np.where(np.iscomplex(unique),np.cos(np.angle(unique)),np.NaN)
    S = np.max(real)<=margin

    # got order to be sorted in
    if sortby == 'F':
        ind = np.argsort(F)
    elif sortby == 'D':
        ind = np.argsort(D)
    else:
        ind = range(len(unique)) 

    # place data in a dict and sort
    res = []
    Mode = 0
    for i in ind:
        res_dict = {}
        res_dict['Real'] = real[i]
        res_dict['Imag'] = imag[i]
        res_dict['Frequency'] = F[i]
        res_dict['Damping'] = D[i]
        res_dict['Stable'] = S
        res_dict['Eigen Vector'] = unique_vecs[i]
        res_dict['Mode'] = Mode
        Mode += 1
        res.append(res_dict)

    return res
