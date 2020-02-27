import sympy as sym
from sympy.abc import t
import numpy as np
from scipy.optimize import minimize
import pandas as pd

# and function to sub in a time dependent variable without altering its derivatives
def msub(expr,v,sub,derivatives):
    """
    Substitutes the symbol 'sub' with the value 'v' in the expression 'expr',
    without changing the first 'derivatives' derivatives of 'sub'.
    i.e. if 'derivatives' = 2 the first two derivatives of 'sub' will be unaffected. 
    All other derivatives of 'sub' will become zero! (hangover of how sympy 
    subing works, as the derivate of a constant is zero)
    """
    temps = list(sym.symbols(f'temps:{derivatives}'))
    derivs = list(sym.symbols(f'temps:{derivatives+1}'))
    derivs[0] = v
    for i in range(1,derivatives+1):
        derivs[i]=derivs[i-1].diff(t)
        
    expr = expr.subs({derivs[i]:temps[i-1] for i in range(1,derivatives+1)})
    expr = expr.subs(v,sub)
    expr = expr.subs({temps[i-1]:derivs[i] for i in range(1,derivatives+1)})
    return expr

def LinearEoM_func(dof2Model,FwtParams,ExtForces,ignore=[]):
    p = FwtParams
    # create complete EoM
    Q = ExtForces.Q
    EoM = dof2Model.X.subs({dof2Model.F[i]:Q[i] for i in range(0,p.qs)})

    # sub all in those to ignore
    tup = p.GetTuple(ignore)
    eqs = (EoM.subs({v:v.value for v in tup}))

    #calculate the jacobian
    jac = eqs.jacobian(p.x)

    # symbols for the fixed point values
    lp = sym.Matrix(sym.symbols(f'lp:{p.qs*2}'))

    # convert displacements
    for i in range(0,p.qs):
        jac = msub(jac,p.q[i],lp[i*2],2)

    # convert velocities
    for i in range(0,p.qs):
        jac = msub(jac,p.qd[i],lp[i*2+1],1)
    
    tup = tuple(ignore)

    # return the function
    return sym.lambdify((*tup,lp),jac)
