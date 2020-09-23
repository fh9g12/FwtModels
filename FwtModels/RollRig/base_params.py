import sys, os
sys.path.insert(1, os.path.join(sys.path[0], '../..'))
import ModelFramework as mf
import sympy as sym
import numpy as np

def base_params(Dofs = 3):
    #2 Dof System
    p = mf.ModelParameters.DynamicModel(Dofs)


    #Geometry Parameters
    p.c = mf.ModelSymbol(value = 0.15,string = 'c') # chord of wing
    p.s = mf.ModelSymbol(value = 1,string = 's')    # span of wing
    p.sigma = mf.ModelSymbol(value = 0.28,string = 'sigma')  # FWT as percentage of total span

    p.s_f = p.s*p.sigma*sym.Rational(1,2) # span of each FWT
    p.s_w = p.s*(1-p.sigma) # span of main wing

    # Mass Parameter
    p.m_w = mf.ModelSymbol(value = 0.22, string = 'm_w') # mass per unit span of inner wing
    p.m_f = mf.ModelSymbol(value = 0.038,string = 'm_f') # mass of each FWT
    p.m_l = mf.ModelSymbol(value = 0.0275,string = 'm_l') # mass of each FWT lock

    p.I_xx_f = mf.ModelSymbol (value = 0.1, string = 'I_xxf')
    p.I_xx_f = sym.Rational(1,12)*p.m_f*p.s_f**2 # inertia of FWT (uniform bar)
    p.I_xx_w = mf.ModelSymbol (value = 0.1, string = 'I_xxw')
    p.I_xx_w = sym.Rational(1,12)*p.m_w*p.s_w**2 # inertia of wing (uniform bar)


    # Attitude Parmas
    p.rho = mf.ModelSymbol(value = 1.225,string = 'rho')                 # density
    p.V = mf.ModelSymbol(value = 10,string = 'V')                        # velocity
    p.g  = mf.ModelSymbol(value = 9.81,string = 'g')                     # gravity
    p.alpha_r = mf.ModelSymbol(value = np.deg2rad(3),string = 'alpha_r') # root AoA
    p.alpha_c = mf.ModelSymbol(value = 0, string = 'alpha_c')            # FWT camber



    p.alpha_0 = mf.ModelSymbol(value = 0,string = 'alpha_0')             # FWT alpha
    p.alphadot_0 = mf.ModelSymbol(value = 0,string = 'alphadot_0')       # FWT alphadot
    p.alpha_1 = mf.ModelSymbol(value = 0,string = 'alpha_1')             # FWT alpha
    p.alphadot_1 = mf.ModelSymbol(value = 0,string = 'alphadot_1')       # FWT alphadot
    p.alpha_2 = mf.ModelSymbol(value = 0,string = 'alpha_2')             # FWT alpha
    p.alphadot_2 = mf.ModelSymbol(value = 0,string = 'alphadot_2')       # FWT alphadot

    p.c_d_max = mf.ModelSymbol(value = 1,string='C_Dmax')
    p.a_0 = mf.ModelSymbol(value = 2*np.pi, string = 'a_0')   # C_L at the root
    p.a_1 = mf.ModelSymbol(value = 2*np.pi, string = 'a_1')   # C_L at the root
    p.beta = mf.ModelSymbol(value = 2*np.pi, string = 'beta' )
    p.eta_0 = mf.ModelSymbol(value = 0,string='eta_0')            # FWT Twist
    p.eta_1 = mf.ModelSymbol(value = 0,string='eta_1')            # FWT Taper

    # fwt params
    p.Lambda = mf.ModelSymbol(value = np.deg2rad(10),string = 'Lambda') # flare angle
    p.y_f = sym.Symbol('y_f')# spanwise location
    p.y_w = sym.Symbol('y_w') # spanwise location

    #Gust Velocity
    p.w_g = mf.ModelSymbol(value = 0,string='w_g')


    ## FWTD Structural Parameters
    p.y_0 = sym.Symbol('y_0')
    p.x_0 = sym.Symbol('x_0')

    p.mode = mf.ModelSymbol(value ='Free',string='Mode')

    ## Numeric Model Constants
    p.fp = mf.ModelMatrix(value =[0]*p.qs*2,symbols=sym.symbols(f'qtilde:{p.qs*2}')) # The stationary point
    return p