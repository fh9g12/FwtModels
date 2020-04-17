import sys, os
sys.path.insert(1, os.path.join(sys.path[0], '../..'))
import ModelFramework as mf
import sympy as sym
import numpy as np

def base_params(generalised_coords):
    p = mf.ModelParameters.DynamicModel(generalised_coords)

    ## Create some global parameters
    p.c = mf.ModelSymbol(value = 1.8,string = 'c') # chord of wing
    p.s_t = mf.ModelSymbol(value = 12,string = 's_t') # total semi-span of wing
    p.rho_t = mf.ModelSymbol(value = 20,string = 'rho_t') # per per unit area (kg/m^2)
    p.ratio_fwt = mf.ModelSymbol(value = 0.2,string = 'sigma') # percentage of wing that is a FWT


    # Respective spans
    p.s_0 =  p.s_t*(1-p.ratio_fwt) # span of inner wing
    p.s_1 =  p.s_t*p.ratio_fwt # span of inner wing

    # Inner Wing stiffness Properties
    p.EI = mf.ModelSymbol(value = 1e7,string = 'EI') # Bending stiffness for the wing
    p.GJ = mf.ModelSymbol(value = 1e6,string = 'GJ') # Torsional Stiffness for the wing

    # Location of flexural axes
    p.e_0 = mf.ModelSymbol(value = 0,string = 'e_0')
    p.x_f0 = sym.Rational(1,4)*p.c + p.e_0*p.c
    p.e_1 = mf.ModelSymbol(value = 0,string = 'e_1') # e for the FWT
    p.x_f1 = sym.Rational(1,4)*p.c + p.e_1*p.c

    # FWT Properties
    p.m_factor = mf.ModelSymbol(value = 1, string = 'delta_m')
    p.m_1 = p.rho_t*p.c*p.s_t*p.ratio_fwt*p.m_factor

    p.I_xx_1 = sym.Rational(1,12)*p.m_1*p.s_1**2 # inertia of FWT (uniform bar)
    p.Lambda = mf.ModelSymbol(value = 0,string = 'Lambda') # Flare Angle
    p.Delta_m = mf.ModelSymbol(value = 0,string = 'Delta_m') # additional mass to apply at the FWT CoM

    # Symbols to translate along inner wing and FWT
    p.y_0 = sym.Symbol('y_0') # inner wing y chord
    p.x_0 = sym.Symbol('x_0') # inner wing x chord
    p.y_1 = sym.Symbol('y_1') # FWT y chord
    p.x_1 = sym.Symbol('x_1') # FWT x chord

    ## Aero Parameters
    p.rho = mf.ModelSymbol(value = 1.225,string = 'rho')                 # density
    p.V = mf.ModelSymbol(value = 10,string = 'V')                        # velocity
    p.g  = mf.ModelSymbol(value = 9.81,string = 'g')                     # gravity
    p.alpha_r = mf.ModelSymbol(value = 0,string = 'alpha_r') # root AoA
    p.M_thetadot = mf.ModelSymbol(value = -1.2,string = 'M_thetadot')    # Unsteady Torsional Term 

    ## Main Wing Specific
    p.a_0 = mf.ModelSymbol(value = 2*np.pi,string = 'a_0')               # C_L slope of main wing

    ## FWT Specific
    p.a_1 = p.a_0 - p.a_0/p.s_1*p.y_1                                    # C_L slope of FWT
    p.alpha_1 = mf.ModelSymbol(value = 0,string = 'alpha_1')             # FWT alpha
    p.alphadot_1 = mf.ModelSymbol(value = 0,string = 'alphadot_1')       # FWT alphadot

    ## Numeric Model Constants
    p.fp = mf.ModelMatrix(value =[0]*p.qs*2,symbols=sym.symbols(f'qtilde:{p.qs*2}')) # The stationary point

    # Factor to change the size of joint values
    p.eta = mf.ModelMatrix(value = [1]*(generalised_coords),symbols=sym.symbols(f'eta:{p.qs}'))
    return p