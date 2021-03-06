import sys, os
sys.path.insert(1, os.path.join(sys.path[0], '../..'))
import moyra as ma
import sympy as sym
import numpy as np

def base_params(Dofs = 3,panels = 20):
    #2 Dof System
    p = ma.DynamicModelParameters(Dofs)


    #Geometry Parameters
    p.c = ma.ModelSymbol(value = 0.15,string = 'c') # chord of wing
    p.s = ma.ModelSymbol(value = 1,string = 's')    # span of wing
    p.sigma = ma.ModelSymbol(value = 0.28,string = 'sigma')  # FWT as percentage of total span

    p.s_f = p.s*p.sigma*sym.Rational(1,2) # span of each FWT
    p.s_w = p.s*(1-p.sigma) # span of main wing

    # Mass Parameter
    p.m_w = ma.ModelSymbol(value = 0.22, string = 'm_w') # mass per unit span of inner wing
    p.m_f = ma.ModelSymbol(value = 0.038,string = 'm_f') # mass of each FWT

    p.I_xxf = ma.ModelSymbol (value = 0.1, string = 'I_xxf')
    p.I_xxw = ma.ModelSymbol (value = 0.1, string = 'I_xxw')


    # Attitude Parmas
    p.rho = ma.ModelSymbol(value = 1.225,string = 'rho')                 # density
    p.V = ma.ModelSymbol(value = 10,string = 'V')                        # velocity
    p.g  = ma.ModelSymbol(value = 9.81,string = 'g')                     # gravity
    p.alpha_r = ma.ModelSymbol(value = np.deg2rad(3),string = 'alpha_r') # root AoA
    p.alpha_c = ma.ModelSymbol(value = 0, string = 'alpha_c')            # FWT camber



    p.alpha_0 = ma.ModelSymbol(value = 0,string = 'alpha_0')             # FWT alpha
    p.alphadot_0 = ma.ModelSymbol(value = 0,string = 'alphadot_0')       # FWT alphadot
    p.alpha_1 = ma.ModelSymbol(value = 0,string = 'alpha_1')             # FWT alpha
    p.alphadot_1 = ma.ModelSymbol(value = 0,string = 'alphadot_1')       # FWT alphadot
    p.alpha_2 = ma.ModelSymbol(value = 0,string = 'alpha_2')             # FWT alpha
    p.alphadot_2 = ma.ModelSymbol(value = 0,string = 'alphadot_2')       # FWT alphadot

    p.c_d_max = ma.ModelSymbol(value = 1,string='C_Dmax')
    p.a_0 = ma.ModelSymbol(value = 2*np.pi, string = 'a_0')   # C_L at the root
    p.a_1 = ma.ModelSymbol(value = 2*np.pi, string = 'a_1')   # C_L at the root
    p.beta = ma.ModelSymbol(value = 2*np.pi, string = 'beta' ) # Aileron Angle

    # fwt params
    p.Lambda = ma.ModelSymbol(value = np.deg2rad(10),string = 'Lambda') # flare angle
    p.y_f = sym.Symbol('y_f')# spanwise location
    p.y_i = sym.Symbol('y_i') # spanwise location

    #Gust Velocity
    p.w_g = ma.ModelSymbol(value = 0,string='w_g')

    #CLosed loop gains for roll control
    p.p = ma.ModelSymbol(value = 1,string='p')

    # Time constnat for torque
    p.T = ma.ModelSymbol(value = 1,string='T')

    ## Numeric Model Constants
    p.a = ma.ModelMatrix(value =[np.pi*2]*panels,length=panels, string='a') # The stationary point

    p.y_w = ma.ModelSymbol(value = 0,string = 'y_w') # y offset of wing CoM
    p.z_w = ma.ModelSymbol(value = 0,string = 'z_w') # z offset of wing CoM
    p.l_f = ma.ModelSymbol(value = 0,string = 'l_f') # moment arm of FWT

    #integration constant
    p.y_n = ma.ModelSymbol(value = 0,string = 'y_n') # integration constant


    ## FWTD Structural Parameters
    p.y_0 = sym.Symbol('y_0')
    p.x_0 = sym.Symbol('x_0')

    ## Numeric Model Constants
    p.fp = ma.ModelMatrix(value =[0]*p.qs*2,length=p.qs*2, string='qtilde') # The stationary point
    return p