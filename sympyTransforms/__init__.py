from .Homogenous import HomogenousTransform,Vee,Wedge

# monkey patch lambdify to use cse
from sympy.utilities.lambdify import _EvaluatorPrinter
from .overload_lambdify import doprint

_EvaluatorPrinter.doprint = doprint