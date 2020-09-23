from .base_params import base_params

import sys, os
sys.path.insert(1, os.path.join(sys.path[0], '../..'))
import ModelFramework as mf
import ModelFramework.ExternalForces as ef

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

def GenRunData_StepTorque(filename,qs,ic,end_time,params,xNames=None,additional_cols = {},sample_freq=100):
    # if no xNames provided create some
    if xNames==None:
        xNames = [f'x{i}' for i in range(qs*2)]
        
    # Set the parameters
    p = base_params(qs)
    p_vars = vars(p)
    for string,value in params.items():
        if string in p_vars:
            p_vars[string].value = value
    # Load the Model   
    sm = mf.SymbolicModel.from_file(filename)
    sm.ExtForces = ef.CompositeForce([sm.ExtForces,ef.CustomForce(None)])
             
    # Create Numeric Model
    nm = mf.NumericModel.from_SymbolicModel(p,sm)
    ext_f = nm.ExtForces.force_funcs[0]
            
    # Create the Torque Force
    torque_period = 0.1
    torque_max = p.beta.value*p.V.value**2

    def torque(tup,x,t):
        if t<3:
            return 0
        elif t<torque_period/2+3:
            return torque_max*0.5*(1- np.cos(2*np.pi*(t-2)/torque_period))
        else:
            return torque_max

    def custom_force(tup,x,t):
        res = np.array([[0.]]*p.qs)
        res[0] = torque(tup,x,t)
        return res
    
    nm.ExtForces.force_funcs = [ext_f , custom_force]
    
    # Solve the problem
    yData = solve_ivp(lambda t,y:nm.deriv(t,y,p.GetNumericTuple(y,t)),(0,end_time),ic,max_step=1/sample_freq)
    
    #interpolate to reduce t and y points
    t = np.linspace(0,end_time,(end_time*sample_freq)+1)    
    int_func = interp1d(yData.t,yData.y)  
    yi = int_func(t)
    
    #generate torque data
    tau = [torque((0,),[0]*p.qs*2,i) for i in t]
    #generate a dictionary of the x data
    q_data = [dict(zip(xNames,row)) for row in yi.T]    

    data = [{'t':_t,'torque':tau[i],**q_data[i],**params,**additional_cols} for i,_t in enumerate(t)]
    return data