from sympy import *
def get_M():
	m_w = Symbol('m_w')
	s = Symbol('s')
	sigma = Symbol('sigma')
	e = MutableDenseMatrix([[m_w*s**2*sigma**2/12 - m_w*s**2*sigma/6 + m_w*s**2/12]])
	return e
def get_f():
	e = ImmutableDenseMatrix([[0]])
	return e
def get_T():
	m_w = Symbol('m_w')
	s = Symbol('s')
	sigma = Symbol('sigma')
	t = Symbol('t')
	q0 = Function('q0')
	e = m_w*s**2*sigma**2*Derivative(q0(t), t)**2/24 - m_w*s**2*sigma*Derivative(q0(t), t)**2/12 + m_w*s**2*Derivative(q0(t), t)**2/24
	return e
def get_U():
	e = 0
	return e
def get_Q():
	V = Symbol('V')
	a_0 = Symbol('a_0')
	c = Symbol('c')
	rho = Symbol('rho')
	s = Symbol('s')
	sigma = Symbol('sigma')
	t = Symbol('t')
	q0 = Function('q0')
	e = ImmutableDenseMatrix([[-V*a_0*c*rho*s**3*(1 - sigma)**3*Derivative(q0(t), t)/24]])
	return e
