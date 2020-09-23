import sympy as sym
import sympy.physics.mechanics as me
from sympy.abc import t as time
import pandas as pd
import seaborn as sns

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import fsolve,least_squares,root

import sys, os
import pickle

sys.path.insert(1, os.path.join(sys.path[0], '../..'))
import custom_plot_objects as cpo

import ModelFramework as mf
import ModelFramework.Elements as ele
import ModelFramework.ExternalForces as ef
import FwtModels.RectWing as rw
import FwtModels.RollRig as RollRig
import multiprocessing as mp
from tqdm import tqdm

from matplotlib.lines import Line2D

def GenRunData(filename,qs,ic,end_time,params,xNames=None,additional_cols = {},sample_freq=100):
    # if no xNames provided create some
    if xNames==None:
        xNames = [f'x{i}' for i in range(qs*2)]
        
    # Set the parameters
    p = RollRig.base_params(qs)
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
        if t<2:
            return 0
        elif t<torque_period/2+2:
            return torque_max*0.5*(1- np.cos(2*np.pi*(t-2)/torque_period))
        else:
            return torque_max

    def custom_force(tup,x,t):
        res = [[0]]*p.qs
        res[0][0] = torque(tup,x,t)
        return np.array(res)
    
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

    data = [{'t':_t,'torque':tau[i],'Locked':True,**q_data[i],**params,**additional_cols} for i,_t in enumerate(t)]
    return data

# Add fixed params to dict
params = {}
params['c'] = 0.067
params['s'] = 1
params['sigma'] = 0.28
params['m_w'] = 0.22
params['m_f'] = 0.038
params['m_l'] = 0#0.0275
params['alpha_r'] = 0
params['c_d_max'] = 1
params['a_0'] = 2*np.pi
params['a_1'] = 2*np.pi
#arams['.tau.value = 0
params['w_g'] = 0
params['beta'] = 0.0036
params['V'] = 22.5
params['Lambda'] = np.deg2rad(10)
params['Mode'] = 'Free'

jobs=[]
# Make jobs for FWT geoms
cambers = np.deg2rad([0,5,10])
twists = np.deg2rad([0,-10,10])
sigmas = [0.28,0.4,0.5]
tapers = [0,-0.5,0.5]

iters = [cambers,twists,sigmas,tapers]
names = ['alpha_c','eta_0','sigma','eta_1']

default = [i[0] for i in iters]
jobs =[]
for i,_iter in enumerate(iters):
    row = default.copy()
    for val in _iter:
        row_0 = row.copy()
        row_0[i] = val
        jobs.append(dict(zip(names,row_0)))

# for each job do multiple velocities
jobs = [{**d,'V':v} for d in jobs for v in np.linspace(15,30,4)] # for each velocity
jobs = [{**d,'beta':v} for d in jobs for v in [0.0036,0.0015]] # for two betas

# make a fixed and free version of each job
jobs = [{**row,**mode}for row in jobs for mode in [{'Mode':'Fixed','m_f':0.038,'a_1':2*np.pi},{'Mode':'Removed','m_f':0,'a_1':0}] ]

print(f'Running {len(jobs)} jobs...')

pool = mp.Pool(mp.cpu_count())
names = ['Roll','Roll Rate','Left FWT Angle','Left FWT Velocity','Right FWT Angle','Right FWT Velocity']

res=[]

pbar = tqdm(total=len(jobs))
def update(*a):
    pbar.update()
    # tqdm.write(str(a))
for job in jobs:
    res.append(pool.apply_async(GenRunData, args = ('/Users/fintan/Documents/GitHub/FwtModels/Workbooks/14_Rolling Rig V2/RollRigModel-Fixed.py',1,[0]*2,10,{**params,**job},names[:2],{n:0 for n in names[2:]},100), callback=update))
pool.close()
pool.join()

ls = []
for i in res:
    ls += i.get()
df = pd.DataFrame(ls)
df.to_pickle('/Users/fintan/Documents/GitHub/FwtModels/Workbooks/14_Rolling Rig V2/FixedData.pkl')

print('Complete!')
