import sympy as sym
import numpy as np
import pandas as pd
from scipy.linalg import eig
from scipy.optimize import fsolve,least_squares
from sympy.physics.mechanics import msubs

import sys,os
sys.path.insert(1, os.path.join(sys.path[0], '../..'))
import ModelFramework as mf

def fixed_point_finder(p,model,vars_ls,jac=True,fixed_point_gen=None,additional_func = {}):
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

    # If caluclating fixed points generate require objective functions
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
    fixed_point_dataframe = []
    qs = []
    for i in range(len(df_perms)):
        values = tuple([v for k,v in df_perms[i]])      
        # Calc fixed point
        #set the initial guess (if v=0 set to FWT dropped doen else use previous result)
        if string_perms[i]["V"] == 0:
            guess = [0]*p.qs
            guess[-1] = np.pi/2
            q = fsolve(lambda q,v: func_obj_v0(q,v)[:,0],guess,factor = 1,args=(values,))  
        else:         
            guess = [0]*p.qs #if len(qs)==0 else qs[-1]
            guess[-1] += 0.1 #if len(qs)==0  else 0 # gets you off the hill (where jacobian is zero for fwt...)

        lb = [-np.inf]*p.qs
        ub = [np.inf]*p.qs
        lb[0] =-0.1
        ub[0] = 0.1
        lb[-1] = -2
        ub[-1] = 2
        bounds = (lb,ub)

        q = least_squares(lambda q,v:func_obj(q,values)[:,0],guess,bounds = bounds,x_scale='jac',
                            method='dogbox',jac=func_jac_obj if jac else '2-point',args = (values,)).x 
        qs.append(q)
        additional = {k:v(q,values) for k,v in additional_func.items()}
        fixed_point_dataframe.append({**string_perms[i],**additional,'q':q})       
    return pd.DataFrame(fixed_point_dataframe)