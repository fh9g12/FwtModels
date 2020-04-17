import FwtModels as fwt
import FwtModels.RectWing as rw
import ModelFramework as mf
import numpy as np


def test_easy():
    b_modes = 3
    t_modes = 3
    model_num = 1
    p = rw.base_params(b_modes + t_modes + 1)
    sm = mf.SymbolicModel.from_file(f'tests/{b_modes}B{t_modes}T-M{model_num}.py')
    
    # Set HALE Specific parameters
    p.rho.value = 0.0889
    p.s_t.value = 16
    p.c.value = 1
    p.EI.value = 2e4
    p.GJ.value = 1e4
    p.rho_t.value = 0.75
    p.e_0.value = 0.25
    p.e_1.value = 0
    
    # set Parameter permutations
    
    vars_ls =[]
    #vars_ls.append((p.m_factor,[0.5,1,1.5]))
    #vars_ls.append((p.Lambda,np.deg2rad([10,17.5,25])))
    #vars_ls.append((p.alpha_r,np.deg2rad([0,5,10])))
    #vars_ls.append((p.ratio_fwt,[0,0.1,0.2,0.3]))
    #vars_ls.append((p.V,np.linspace(0,200,201))) # ensure velocity last so that fixed points iterats up the velocity

    vars_ls.append((p.m_factor,[1]))
    vars_ls.append((p.Lambda,np.deg2rad([10])))
    vars_ls.append((p.alpha_r,np.deg2rad([0,5,10])))
    vars_ls.append((p.ratio_fwt,[0,0.1,0.2,0.3]))
    vars_ls.append((p.V,np.linspace(0,40,80))) # ensure velocity last so that fixed points iterats up the velocity
    
    # generate result
    res = rw.eigen_perm_params(p,sm,vars_ls,np.isin(model_num,[1,2,3,4]))
    assert (3==3)