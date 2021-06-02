from .base_params import base_params

import sys, os
sys.path.insert(1, os.path.join(sys.path[0], '../..'))
import moyra as ma
import moyra.forces as ef

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from scipy.optimize import minimize


def calc_coast_angle(filename, qs, ic, params, xNames=None, additional_cols={}):
    # if no xNames provided create some
    if xNames is None:
        xNames = [f'x{i}' for i in range(qs*2)]
        
    # Set the parameters
    p = base_params(qs)
    p_vars = vars(p)
    for string, value in params.items():
        if string in p_vars:
            p_vars[string].value = value
    # Load the Model   
    sm = ma.SymbolicModel.from_file(filename)
    sm.ExtForces = ef.CompositeForce([sm.ExtForces, ef.Customaorce(None)])
             
    # Create Numeric Model
    nm = ma.NumericModel.from_SymbolicModel(p,sm)
    ext_f = nm.ExtForces.force_funcs[0]

    # Create Objetive Function
    def objective_func(x, roll, roll_rate, index):
        y=[0]*6
        y[0] = roll
        y[1] = roll_rate
        y[index-1] = np.deg2rad(x[0])
        res = nm.deriv(y, p.GetNumericTuple(y, 0),0)
        return res[index]**2
    
    # create function to find a good inital guess
    folds = np.linspace(-60, 60, 15)
    f = lambda x,index : [objective_func([x_i], ic[0], ic[1], index) for x_i in x]

    # find left and right cruise angle
    res_left = minimize(objective_func, [folds[np.argmin(f(folds, 3))]], args=(ic[0], ic[1], 3), bounds=((-90, 90),))
    res_right = minimize(objective_func, [folds[np.argmin(f(folds, 5))]], args=(ic[0], ic[1], 5), bounds=((-90, 90),))

    # calculate roll torque in coast postion
    if res_left.success and res_right.success:
        y = [ic[0],ic[1], np.deg2rad(res_left.x[0]), 0, np.deg2rad(res_right.x[0]), 0]
        torq = nm.deriv(y, p.GetNumericTuple(y, 0),0)[1]
    else:
        torq = np.nan

    # return result
    result = []
    result.append({"Torque":torq,"Side":"Left",
                "CoastAngle":res_left.x[0] if res_left.success else np.nan,
                **params, **additional_cols })
    result.append({"Torque":torq,"Side":"Right",
                "CoastAngle":res_right.x[0] if res_right.success else np.nan,
                **params, **additional_cols})
    return result

def Calc_coast(numeric_model,p,ic):
    # Create Objetive Function
    def objective_func(x, ic):
        y=[0]*6
        y[0] = ic[0]
        y[1] = ic[1]
        y[2] = x[0]
        y[4] = x[1]
        tup = p.GetNumericTuple(y, 0)  
        forces = -numeric_model.f(y,tup)+numeric_model.ExtForces(y,tup,0)
        return forces[1][0]**2 + forces[2][0]**2
    # find left and right cruise angle
    res = minimize(objective_func, [0,0], args=(ic,))
    if res.success:
        return [ic[0],ic[1],res.x[0],0,res.x[1],0]
    raise Exception('Failed to calulate coast angles')

def GenRunData_StepTorque(filename, qs, ic, end_time, params,panels=20, calc_coast=False, xNames=None, additional_cols={}, sample_freq=100,events = None):
    # if no xNames provided create some
    if xNames is None:
        xNames = [f'x{i}' for i in range(qs*2)]
        
    # Set the parameters
    (sm,p) = ma.SymbolicModel.from_file(filename)
    p_vars = vars(p)
    for string, value in params.items():
        if string in p_vars:
            p_vars[string].value = value
    # Load the Model 
    
    sm.ExtForces = ef.CompositeForce([sm.ExtForces, ef.CustomForce(None)])
             
    # Create Numeric Model
    nm = ma.NumericModel.from_SymbolicModel(p, sm)
    ext_f = nm.ExtForces.force_funcs[0]
    
    # calcualte coast angles
    ic = Calc_coast(nm,p,ic) if calc_coast else ic
    print(ic)

    # Create the Torque Force
    torque_period = p.T.value
    torque_max = -p.beta.value

    def torque(x, tup, t):
        res = np.array([[0.]]*int(len(x)/2))
        if t<torque_period:
            res[0] = t/torque_period*torque_max
        else:
            res[0] = torque_max          
        return res

    def fold_limit(x, tup, t):
        
        res = np.array([[0.]]*int(len(x)/2))
        lim = np.deg2rad(120)
        res[1] = 0 if abs(x[2]) < lim else -(x[2]%lim)*0.3
        res[2] = 0 if abs(x[4]) < lim else (x[4]%lim)*0.3
        return res
    
    if qs == 1:
        nm.ExtForces.force_funcs = [ext_f, torque]
    else:
        nm.ExtForces.force_funcs = [ext_f, lambda x,tup,t: torque(x,tup,t)+fold_limit(x,tup,t)]
    
    # Solve the problem
    t = np.linspace(0, end_time, (end_time*sample_freq)+1)
    y_data = solve_ivp(lambda t, y: nm.deriv(y,p.GetNumericTuple(y, t),t), (0, end_time), ic,t_eval = t,events = events)
    
    #interpolate to reduce t and y points
    
    t = y_data.t
    yi = y_data.y

    #generate torque data
    tau = [torque([0]*p.qs*2,None, i)[0][0] for i in t]

    #generate a dictionary of the x data
    q_data = [dict(zip(xNames, row)) for row in yi.T]
    data = [{'t':_t, 'torque':tau[i],#'lift_mom':hinge_moment_lift[i],'right_lift':right_lift[i],'left_lift':left_lift[i],
            #'centri_mom':hinge_moment_centri[i],'grav_mom':hinge_moment_grav[i],
            **q_data[i], **params, **additional_cols} for i, _t in enumerate(t)]
    return data