import sys, os
sys.path.insert(1, os.path.join(sys.path[0], '../..'))
import moyra as ma
import sympy as sym
import numpy as np

def base_params(panels = 10,dofs=2):
    #2 Dof System
    p = ma.DynamicModelParameters(dofs)

    ## Create some global parameters
    p.c = ma.ModelSymbol(value = 0.15,string = 'c') # chord of wing

    p.m = ma.ModelSymbol(value = 1,string = 'm')      # mass of fwt
    p.m_1 = ma.ModelSymbol(value = 1, string='m_1')   # mass of additional weight
    p.m_w = ma.ModelSymbol(value = 4,string = 'm_w')  # mass of inner wing
    p.s = ma.ModelSymbol(value = 1,string = 's')      # span of fwt

    p.f_0 = ma.ModelSymbol(value = 2,string = 'f_0')  # the frequency of the first bending mode 
    p.k_w = ma.ModelSymbol(value = 2,string = 'k_w')  # Stiffness of the inner wing
    p.c_w = ma.ModelSymbol(value = 2,string = 'c_w')  # Inner wing rotation coefficent
    p.d_w = ma.ModelSymbol(value = 2, string='d_w') # Damping coefficent for inner wing
    p.d_a = ma.ModelSymbol(value = 2, string='d_a') # Aero Damping coefficent for inner wing
    p.k_fwt = ma.ModelSymbol(value = 0,string = 'k_fwt')  # Stiffness of the folding wing tip

    p.I_xx = ma.ModelSymbol(value = 1,string = 'I_xx')

    p.l_com = ma.ModelSymbol(value = 0.5,string='l_com') # distance from Hinge to CoM
    p.l_m = ma.ModelSymbol(value = 0.5,string='l_m') # distance from Hinge to extra weight

    p.Lambda = ma.ModelSymbol(value = np.deg2rad(10),string = 'Lambda') # flare angle

    p.rho = ma.ModelSymbol(value = 1.225,string = 'rho')                 # Density
    p.V = ma.ModelSymbol(value = 10,string = 'V')                        # Velocity
    p.g  = ma.ModelSymbol(value = 9.81,string = 'g')                     # Gravity
    p.alpha_r = ma.ModelSymbol(value = np.deg2rad(3),string = 'alpha_r') # Root AoA

    p.M_thetadot = ma.ModelSymbol(value = 1.2,string = 'M_thetadot')     # Unsteady Term
    p.alpha_1 = ma.ModelSymbol(value = 0,string = 'alpha_1')             # FWT alpha
    p.alphadot_1 = ma.ModelSymbol(value = 0,string = 'alphadot_1')       # FWT alphadot
    p.clip_factor = ma.ModelSymbol(value = 100, string = 'mu') # roundedness of C_l Curve
    p.c_d_max = ma.ModelSymbol(value = 1,string='C_Dmax')


    #Gust Velocity
    p.w_g = ma.ModelSymbol(value = 0,string='w_g')

    ## FWTD Structural Parameters
    p.y_0 = sym.Symbol('y_0')
    p.y_i = sym.Symbol('y_i')
    p.x_0 = sym.Symbol('x_0')

    ## Numeric Model Constants
    p.fp = ma.ModelMatrix(value =[0]*p.qs*2,length=p.qs*2, string='qtilde') # The stationary point
    p.a = ma.ModelMatrix(value =[np.pi*2]*panels,length=panels, string='a') # The aero derivative on each panel

    ## FWT C_L 
    p.a_0 = ma.ModelSymbol(value = 2*np.pi,string = 'a_0')   # C_L at the root
    return p