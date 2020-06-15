import sys, os
sys.path.insert(1, os.path.join(sys.path[0], '../..'))
import ModelFramework as mf
import sympy as sym
import numpy as np

def base_params():
    #2 Dof System
    p = mf.ModelParameters.DynamicModel(2)

    ## Create some global parameters
    p.c_root = mf.ModelSymbol(value = 0.15,string = 'c_root') # chord of wing
    p.c_tip = mf.ModelSymbol(value = 0.15,string = 'c_root') # chord of wing

    p.m = mf.ModelSymbol(value = 1,string = 'm')      # mass of fwt
    p.m_w = mf.ModelSymbol(value = 4,string = 'm_w')  # mass of inner wing
    p.s = mf.ModelSymbol(value = 1,string = 's')      # chord of fwt

    p.f_0 = mf.ModelSymbol(value = 2,string = 'f_0')  # the frequency of the first bending mode 
    p.k_w = mf.ModelSymbol(value = 2,string = 'k_w')  # Stiffness of the inner wing
    p.k_fwt = mf.ModelSymbol(value = 0,string = 'k_fwt')  # Stiffness of the folding wing tip
    #p.k_w = (p.f_0*2*sym.pi)**2*(p.m_w+p.m)             # Stiffness of the inner wing

    p.I_xx = mf.ModelSymbol(value = 1,string = 'I_xx')
    #p.I_xx = sym.Rational(1,12)*p.m*p.s**2*1          # FWT polar moment of inetia
    #p.I_xx = sym.Rational(6.5,81)*p.m*p.s**2 

    p.l_com = mf.ModelSymbol(value = 0.5,string='l_com')

    p.Lambda = mf.ModelSymbol(value = np.deg2rad(10),string = 'Lambda') # flare angle

    p.rho = mf.ModelSymbol(value = 1.225,string = 'rho')                 # density
    p.V = mf.ModelSymbol(value = 10,string = 'V')                        # velocity
    p.g  = mf.ModelSymbol(value = 9.81,string = 'g')                     # gravity
    p.alpha_r = mf.ModelSymbol(value = np.deg2rad(3),string = 'alpha_r') # root AoA
    p.M_thetadot = mf.ModelSymbol(value = 1.2,string = 'M_thetadot')     # Unsteady Term

    p.alpha_1 = mf.ModelSymbol(value = 0,string = 'alpha_1')             # FWT alpha
    p.alphadot_1 = mf.ModelSymbol(value = 0,string = 'alphadot_1')       # FWT alphadot
    p.clip_factor = mf.ModelSymbol(value = 100, string = 'mu') # roundedness of C_l Curve
    p.c_d_max = mf.ModelSymbol(value = 1,string='C_Dmax')


    #Gust Velocity
    p.w_g = mf.ModelSymbol(value = 0,string='w_g')


    ## FWTD Structural Parameters
    p.y_0 = sym.Symbol('y_0')
    p.x_0 = sym.Symbol('x_0')


    ## Numeric Model Constants
    p.fp = mf.ModelMatrix(value =[0]*p.qs*2,symbols=sym.symbols(f'qtilde:{p.qs*2}')) # The stationary point

    ## FWT C_L 
    p.a_0 = mf.ModelSymbol(value = 2*np.pi,string = 'a_0')   # C_L at the root
    p.a = p.a_0 - p.a_0/p.s*p.y_0                            # C_L distrobution
    return p