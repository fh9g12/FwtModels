import pandas as pd
import numpy as np
import sys, os
sys.path.insert(1, os.path.join(sys.path[0], '..'))
import ModelFramework as mf
import FwtModels.RectWing as rw



def JEC2():
    p = rw.base_params(b_modes + t_modes + 1)
    p.rho.value = 1.225
    p.s_t.value = 1.345
    p.c.value = 0.15
    p.EI.value = 61
    p.GJ.value = 200
    p.rho_t.value = 10
    p.e_0.value = 0
    p.e_1.value = 0
    p.ratio_fwt.value = 0.2
    p.alpha_s.value = np.deg2rad(15)
    return p

def JEC():
    p = rw.base_params(b_modes + t_modes + 1)
    p.rho.value = 1.225
    p.s_t.value = 12
    p.c.value = 1.8
    p.EI.value = 9.77e6
    p.GJ.value = 0.99e6
    p.rho_t.value = 19.53
    p.e_0.value = 0.08
    p.e_1.value = 0   
    p.ratio_fwt.value = 0.2
    p.alpha_s.value = np.deg2rad(15)
    return p

def HALE():
    p = rw.base_params(b_modes + t_modes + 1)
    p.rho.value = 0.0889
    p.s_t.value = 16
    p.c.value = 1
    p.EI.value = 2e4
    p.GJ.value = 1e4
    p.rho_t.value = 0.75
    p.e_0.value = 0.25
    p.e_1.value = 0
    p.ratio_fwt.value = 0.2
    p.alpha_s.value = np.deg2rad(15)
    return p




def Model_Eigen(model_num,b_modes,t_modes):
    print(f'Genrating data for model {model_num}')
    try:
        p = JEC()
        dataset_name = 'JEC'
        
        sm = mf.SymbolicModel.from_file(f'{b_modes}B{t_modes}T-M{model_num}.py')    
        vars_ls =[]
        vars_ls.append((p.Lambda,np.deg2rad([10,17.5,25])))
        #vars_ls.append((p.V,np.linspace(0,40,81))) # V must be second
        vars_ls.append((p.V,np.linspace(0,150,151))) # V must be second
        vars_ls.append((p.alpha_r,np.deg2rad([0,5,10])))
        vars_ls.append((p.c_dmax,[0,0.5,1,1.5]))
        #vars_ls.append((p.ratio_fwt,[0,0.1,0.2,0.3]))
        #vars_ls.append((p.ratio_DL,[0,0.05,0.1,0.2]))
        vars_ls.append((p.m_factor,[0.5,1,1.5]))
        
        calc_fixed = True if np.isin(model_num,np.array([1,2,3,4,5])) else False
        flutdf = rw.eigen_perm_params(p,sm,vars_ls,calc_fixed,jac = False)   
        flutdf.to_pickle(f'Eigen_{b_modes}B{t_modes}T-M{model_num}_{dataset_name}.pkl')
        print(f'Genrated data for model {model_num}')
    except:
        print(f'Model {model_num} exited with an error')

b_modes = 3
t_modes = 3
import multiprocessing as mp

pool = mp.Pool(mp.cpu_count())

for k in range(8):
#for k in [5,7]:
    pool.apply_async(Model_Eigen,args=(k,b_modes,t_modes))
pool.close()
pool.join()