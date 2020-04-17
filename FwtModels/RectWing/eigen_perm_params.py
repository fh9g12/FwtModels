import sympy as sym
import numpy as np
import pandas as pd
from scipy.linalg import eig
from scipy.optimize import fsolve
from sympy.physics.mechanics import msubs

import sys,os
sys.path.insert(1, os.path.join(sys.path[0], '../..'))
import ModelFramework as mf

def eigen_perm_params(p,model,vars_ls,calc_fixed_points):
    """
    Method to generate the flutter results for a model for each permutation of the parms in param_perms:
    p - instance of Model Parameters
    model - instance of symbolic model using params defined in p
    param_perms - a list of tuples, each tuple consists of (Symbol, list of values)
    calc_fixed_points - if True, will calc the fixed point for each permutation
    """

    #get list of symbols to key in model
    variables = [k for k,v in vars_ls]

    model_mini = model.msubs(p.GetSubs(0,p.fp,ignore=variables))
    model_lin = model_mini.linearise(p)
    fp_v = p.fp[1::2]
    model_mini_lin = model_lin.msubs(p.GetSubs(0,p.fp,ignore=variables)).msubs({i:0 for i in fp_v})

    # get eigen Matrices and turn into a function
    K,M = model_mini_lin.GeneralEigenProblem(p)
    func = sym.lambdify((variables+[p.fp[::2]]),(K,M))

    def func_no_nan(*args):
        k,m = func(*args)
        return np.where(np.isnan(k),0,k),np.where(np.isnan(m),0,m)

    # If caluclating fixed points generate require objective functions
    if calc_fixed_points:
        f = msubs((model_mini.f-model_mini.ExtForces.Q()),{i:0 for i in p.qd})
        func_obj = sym.lambdify((p.q,variables),f)
        func_jac_obj = sym.lambdify((p.q,variables),f.jacobian(p.q),"numpy")


    # Get all possible combinations of the variables
    perms = np.array(np.meshgrid(*[v for k,v in vars_ls ])).T.reshape(-1,len(vars_ls))

    df_perms = [[(vars_ls[i][0],row[i]) for i in range(len(row))] for row in perms]
    string_perms = [{vars_ls[i][0].name:row[i] for i in range(len(row))} for row in perms]

    #Calc freqs and dampings
    flutdfv2 = []
    qs = []
    for i in range(len(df_perms)):      
        # Calc fixed point
        #set the initial guess (if v=0 set to FWT dropped doen else use previous result)
        if calc_fixed_points:
            if string_perms[i]["V"] == 0:
                guess = [0]*p.qs
                guess[-1] = np.pi/2
            elif i == 0:
                guess = [0]*p.qs
            else:
                guess = qs[-1]          
            values = tuple([v for k,v in df_perms[i]])
            q = fsolve(lambda q,v: func_obj(q,values)[:,0],guess,fprime = func_jac_obj ,factor = 1,args=(values,))    
        else:
            q=[0]*p.qs

        qs.append(q)
        evals, evecs = eig(*func_no_nan(*values,q))

        jac_dat = mf.ExtractEigenValueData_list(evals,evecs,sortby='F')
        #create q and perm data to match size
        for j in range(len(jac_dat)):    
            flutdfv2.append({**jac_dat[j],**string_perms[i],'q':q})
            
    return pd.DataFrame(flutdfv2)

