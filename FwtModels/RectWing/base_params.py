import sys, os
sys.path.insert(1, os.path.join(sys.path[0], '../..'))
import moyra as ma
import sympy as sym
import numpy as np

def base_params(generalised_coords):
    p = ma.DynamicModelParameters(generalised_coords)
    
    ## Create some global parameters
    p.c = ma.ModelSymbol(value = 1.8,string = 'c') # chord of wing
    p.s_t = ma.ModelSymbol(value = 12,string = 's_t') # total semi-span of wing
    p.rho_t = ma.ModelSymbol(value = 20,string = 'rho_t') # per per unit area (kg/m^2)
    p.ratio_fwt = ma.ModelSymbol(value = 0.2,string = 'sigma') # percentage of wing that is a FWT

    p.c_dmax = ma.ModelSymbol(value = 0,string = 'c_dmax') # percentage of wing that is a FWT
    p.clip_factor = ma.ModelSymbol(value = 100, string = 'mu') # roundedness of C_l Curve
    p.alpha_s = ma.ModelSymbol(value = np.deg2rad(15),string = 'alpha_s') # stall Angle


    # Respective spans
    p.s_0 =  p.s_t*(1-p.ratio_fwt) # span of inner wing
    p.s_1 =  p.s_t*p.ratio_fwt # span of inner wing

    # Inner Wing stiffness Properties
    p.EI = ma.ModelSymbol(value = 1e7,string = 'EI') # Bending stiffness for the wing
    p.GJ = ma.ModelSymbol(value = 1e6,string = 'GJ') # Torsional Stiffness for the wing

    # Location of flexural axes
    p.e_0 = ma.ModelSymbol(value = 0,string = 'e_0')
    p.x_f0 = sym.Rational(1,4)*p.c + p.e_0*p.c
    p.e_1 = ma.ModelSymbol(value = 0,string = 'e_1') # e for the FWT
    p.x_f1 = sym.Rational(1,4)*p.c + p.e_1*p.c

    # FWT Properties
    p.m_factor = ma.ModelSymbol(value = 1, string = 'delta_m')
    p.m_1 = p.rho_t*p.c*p.s_t*p.ratio_fwt*p.m_factor

    p.I_xx_1 = sym.Rational(1,12)*p.m_1*p.s_1**2 # inertia of FWT (uniform bar)
    p.Lambda = ma.ModelSymbol(value = 0,string = 'Lambda') # Flare Angle
    p.Delta_m = ma.ModelSymbol(value = 0,string = 'Delta_m') # additional mass to apply at the FWT CoM

    # Symbols to translate along inner wing and FWT
    p.y_0 = sym.Symbol('y_0') # inner wing y chord
    p.x_0 = sym.Symbol('x_0') # inner wing x chord
    p.y_1 = sym.Symbol('y_1') # FWT y chord
    p.x_1 = sym.Symbol('x_1') # FWT x chord

    ## Aero Parameters
    p.rho = ma.ModelSymbol(value = 1.225,string = 'rho')                 # density
    p.V = ma.ModelSymbol(value = 10,string = 'V')                        # velocity
    p.g  = ma.ModelSymbol(value = 9.81,string = 'g')                     # gravity
    p.alpha_r = ma.ModelSymbol(value = 0,string = 'alpha_r') # root AoA
    p.M_thetadot = ma.ModelSymbol(value = -1.2,string = 'M_thetadot')    # Unsteady Torsional Term 

    ## Main Wing Specific
    p.a_0 = ma.ModelSymbol(value = 2*np.pi,string = 'a_0')               # C_L slope of main wing

    ## FWT Specific
    p.a_1 = p.a_0 - p.a_0/p.s_1*p.y_1                                    # C_L slope of FWT
    p.alpha_1 = ma.ModelSymbol(value = 0,string = 'alpha_1')             # FWT alpha
    p.alphadot_1 = ma.ModelSymbol(value = 0,string = 'alphadot_1')       # FWT alphadot

    ## Numeric Model Constants
    p.fp = ma.ModelMatrix(value =[0]*p.qs*2,length=p.qs*2, string='qtilde') # The stationary point

    # Factor to change the size of joint values
    p.eta = ma.ModelMatrix(value = [1]*(generalised_coords),length=generalised_coords, string='eta')

    #Gust Velocity
    p.w_g = ma.ModelSymbol(value = 0,string='w_g')

    p.tau_1 = ma.ModelSymbol(value = 0,string="tau_1")# rate of geometric twist in the FWT

    #FWT Stiffness
    p.k_fwt = ma.ModelSymbol(value = 0, string = 'k_fwt')

    # Aircraft Mass and spring constant
    p.m_ac = ma.ModelSymbol(value = 100,string = 'm_ac')
    p.k_ac = ma.ModelSymbol(value = 1e6, string = 'k_ac')

    # Aircraft Yaw Angle
    p.yaw = ma.ModelSymbol(value = 0, string = 'theta')

    # Locked coast angle
    p.theta_c = ma.ModelSymbol(value = 0, string = 'theta_c')

    return p

def JEC_params(generalised_coords):
    p = base_params(generalised_coords)
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
    p.eta.value = [1]*p.qs
    p.m_ac.value = 1e6
    p.k_ac.value = 1e6
    p.k_fwt.value = 0
    return p

def WT_params(generalised_coords):
    p = base_params(generalised_coords)
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
    p.eta.value = [1]*p.qs
    p.m_ac.value = 1e6
    p.k_ac.value = 1e6
    return p

def HALE_params(generalised_coords):
    p = base_params(generalised_coords)
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
    p.eta.value = [1]*p.qs
    p.m_ac.value = 1e6
    p.k_ac.value = 1e6
    return p