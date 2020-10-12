from .base_params import base_params

import sys, os
sys.path.insert(1, os.path.join(sys.path[0], '../..'))
import ModelFramework as mf
import ModelFramework.ExternalForces as ef

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
    sm = mf.SymbolicModel.from_file(filename)
    sm.ExtForces = ef.CompositeForce([sm.ExtForces, ef.CustomForce(None)])
             
    # Create Numeric Model
    nm = mf.NumericModel.from_SymbolicModel(p,sm)
    ext_f = nm.ExtForces.force_funcs[0]

    # Create Objetive Function
    def objective_func(x, roll, roll_rate, index):
        y=[0]*6
        y[0] = roll
        y[1] = roll_rate
        y[index-1] = np.deg2rad(x[0])
        res = nm.deriv(0, y, p.GetNumericTuple(y, 0))
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
        torq = nm.deriv(0, y, p.GetNumericTuple(y, 0))[1]
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




def GenRunData_StepTorque(filename, qs, ic, end_time, params, xNames=None, additional_cols={}, sample_freq=100):
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
    sm = mf.SymbolicModel.from_file(filename)
    sm.ExtForces = ef.CompositeForce([sm.ExtForces, ef.CustomForce(None)])
             
    # Create Numeric Model
    nm = mf.NumericModel.from_SymbolicModel(p, sm)
    ext_f = nm.ExtForces.force_funcs[0]
            
    # Create the Torque Force
    torque_period = 0.1
    torque_max = p.beta.value

    def torque(tup, x, t):
        if t < 3:
            return -x[0]*30
        elif t < torque_period/2+3:
            return torque_max*0.5*(1- np.cos(2*np.pi*(t-2)/torque_period))
        else:
            return torque_max

    def fold_limit(tup, x, t):
        res = np.array([[0.]]*p.qs)
        lim = np.deg2rad(120)
        res[1] = 0 if abs(x[2]) < lim else -(x[2]%lim)*0.03
        res[2] = 0 if abs(x[4]) < lim else (x[4]%lim)*0.03
        return res


    def custom_force(tup, x, t):
        res = np.array([[0.]]*p.qs)
        res[0] = torque(tup, x, t)
        res += fold_limit(tup, x, t)
        return res
    
    nm.ExtForces.force_funcs = [ext_f, custom_force]
    
    # Solve the problem
    y_data = solve_ivp(lambda t, y: nm.deriv(t, y,p.GetNumericTuple(y, t)), (0, end_time), ic, max_step=1/sample_freq)
    
    #interpolate to reduce t and y points
    t = np.linspace(0, end_time, (end_time*sample_freq)+1)
    int_func = interp1d(y_data.t, y_data.y)
    yi = int_func(t)

    #generate torque data
    tau = [torque((0,), [0]*p.qs*2, i) for i in t]
    #generate a dictionary of the x data
    q_data = [dict(zip(xNames, row)) for row in yi.T]
    data = [{'t':_t, 'torque':tau[i], **q_data[i], **params, **additional_cols} for i, _t in enumerate(t)]
    return data