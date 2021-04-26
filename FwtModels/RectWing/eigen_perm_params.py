import sympy as sym
import numpy as np
import pandas as pd
from scipy.linalg import eig
from scipy.optimize import fsolve,least_squares,root
from sympy.physics.mechanics import msubs

import sys,os
sys.path.insert(1, os.path.join(sys.path[0], '../..'))

def eigen_perm_params(p,model,vars_ls,calc_fixed_points,jac=True,fixed_point_gen=None,fp=None,sortby='F'):
    """
    Method to generate the flutter results for a model for each permutation of the parms in param_perms:
    p - instance of Model Parameters
    model - instance of symbolic model using params defined in p
    param_perms - a list of tuples, each tuple consists of (Symbol, list of values)
    calc_fixed_points - if True, will calc the fixed point for each permutation
    """
    if fp is None:
        fp = [0]*p.qs

    #get list of symbols to key in model
    variables = [k for k,v in vars_ls]

    model_mini = model.msubs(p.GetSubs(0,p.fp,ignore=variables))

    # get eigen Matrices and turn into a function
    K,M = model_mini.GeneralEigenProblemLin(p)

    #get free body problem
    K_v0 = msubs(K,{sym.Symbol(p.V.name):0})
    M_v0 = msubs(M,{sym.Symbol(p.V.name):0})
    func = sym.lambdify((variables,p.fp),(K,M))
    func_v0 = sym.lambdify((variables,p.fp),(K_v0,M_v0))

    # If caluclating fixed points generate require objective functions
    if calc_fixed_points:
        if fixed_point_gen is None:
            f = msubs((model_mini.f-model_mini.ExtForces.Q()),{i:0 for i in p.qd})
        else:
            f = fixed_point_gen(model_mini)
        f_v0 = msubs(f,{sym.Symbol(p.V.name):0})

        func_obj = sym.lambdify((p.q,variables),f,"numpy")
        func_obj_v0 = sym.lambdify((p.q,variables),f_v0,"numpy")

        if jac:
            func_jac_obj = sym.lambdify((p.q,variables),f.jacobian(p.q),"numpy")

    # Get all possible combinations of the variables
    perms = np.array(np.meshgrid(*[v for k,v in vars_ls ])).T.reshape(-1,len(vars_ls))

    df_perms = [[(vars_ls[i][0],row[i]) for i in range(len(row))] for row in perms]
    string_perms = [{vars_ls[i][0].name:row[i] for i in range(len(row))} for row in perms]

    #Calc freqs and dampings
    flutdfv2 = []
    qs = []
    for i in range(len(df_perms)):
        values = tuple([v for k,v in df_perms[i]])      
        # Calc fixed point
        #set the initial guess (if v=0 set to FWT dropped doen else use previous result)
        if calc_fixed_points:
            if string_perms[i]["V"] == 0:
                guess = [0]*p.qs
                guess[-1] = -np.pi/2
                q = fsolve(lambda q,v: func_obj_v0(q,v)[:,0],guess,factor = 1,args=(values,))  
            else:
                if i>0:
                    guess = qs[i-1]
                else:
                    guess = [0]*p.qs
                    guess[-1] =0.1
                q = root(lambda q,v:func_obj(q,values)[:,0],guess,jac=func_jac_obj if jac else None,args = (values,)).x         
        else:
            q=fp

        qs.append(q)
        x = [0]*p.qs*2
        x[::2] = q
        if string_perms[i]["V"] == 0:
            evals, evecs = eig(*func_v0(values,x))
        else:
            evals, evecs = eig(*func(values,x))

        jac_dat = ma.ExtractEigenValueData_list(evals,evecs,sortby=sortby)
        #create q and perm data to match size
        for j in range(len(jac_dat)):    
            flutdfv2.append({**jac_dat[j],**string_perms[i],'q':q})
            
    return pd.DataFrame(flutdfv2)

