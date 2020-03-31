from sympy.utilities.lambdify import _EvaluatorPrinter
from sympy.core.compatibility import (exec_, is_sequence, iterable,
    NotIterable, builtins)
    
def doprint(self, funcname, args, expr):
        """Returns the function definition code as a string."""
        from sympy import Dummy,cse,Symbol

        funcbody = []

        if not iterable(args):
            args = [args]

        argstrs, expr = self._preprocess(args, expr)

        ## --------------- Addition -----------------
        replacments, exprs = cse(expr,symbols=(Symbol(f'rep_{i}')for i in range(1000)))
        if isinstance(expr,tuple):
            expr = tuple(exprs)
        elif isinstance(expr,list):
            expr = exprs
        else:
            expr = exprs[0]
        ## --------------- Addition -----------------

        # Generate argument unpacking and final argument list
        funcargs = []
        unpackings = []

        for argstr in argstrs:
            if iterable(argstr):
                funcargs.append(self._argrepr(Dummy()))
                unpackings.extend(self._print_unpacking(argstr, funcargs[-1]))
            else:
                funcargs.append(argstr)
        arg_vars = ', '.join(funcargs)
        funcsig = f'def {funcname}({arg_vars}):'

        # Wrap input arguments before unpacking
        funcbody.extend(self._print_funcargwrapping(funcargs))

        funcbody.extend(unpackings)

        ## --------------- Addition -----------------
        for variable, expression in replacments:
            funcbody.append(f'{variable} = {self._exprrepr(expression)}')
        ## --------------- Addition -----------------

        funcbody.append('return ({})'.format(self._exprrepr(expr)))

        funclines = [funcsig]
        funclines.extend('    ' + line for line in funcbody)

        return '\n'.join(funclines) + '\n'

def msub(expr,v,sub,derivatives):
    """
    Substitutes the symbol 'sub' with the value 'v' in the expression 'expr',
    without changing the first 'derivatives' time derivatives of 'sub'.
    i.e. if 'derivatives' = 2 the first two time derivatives of 'sub' will be unaffected. 
    All other derivatives of 'sub' will become zero! (hangover of how sympy 
    subing works, as the derivate of a constant is zero)
    """
    from sympy import symbols
    from sympy.abc import t 
    temps = list(symbols(f'temps:{derivatives}'))

    # get a symbol representing each symbol we need to 'save'
    derivs = []
    derivs.append(sub)
    for i in range(derivatives):
        derivs.append(derivs[-1].diff(t))
    derivs = derivs[1:]

    # replace the symbols with temparay symbols
    expr = expr.subs({derivs[i]:temps[i] for i in range(derivatives)})

    # make the actual substitiution
    expr = expr.subs(sub,v)

    # sub back in the upper derivatives
    expr = expr.subs({temps[i]:derivs[i] for i in range(derivatives)})
    return expr