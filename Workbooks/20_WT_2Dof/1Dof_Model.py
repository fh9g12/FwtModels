from sympy import *
def get_M():
	l_com = Symbol('l_com')
	m = Symbol('m')
	l_m = Symbol('l_m')
	m_1 = Symbol('m_1')
	e = MutableDenseMatrix([[l_com**2*m + l_m**2*m_1]])
	return e
def get_f():
	g = Symbol('g')
	l_com = Symbol('l_com')
	m = Symbol('m')
	t = Symbol('t')
	l_m = Symbol('l_m')
	m_1 = Symbol('m_1')
	q0 = Function('q0')
	e = ImmutableDenseMatrix([[g*l_com*m*cos(q0(t)) + g*l_m*m_1*cos(q0(t))]])
	return e
def get_T():
	l_com = Symbol('l_com')
	m = Symbol('m')
	t = Symbol('t')
	l_m = Symbol('l_m')
	m_1 = Symbol('m_1')
	q0 = Function('q0')
	e = l_com**2*m*Derivative(q0(t), t)**2/2 + l_m**2*m_1*Derivative(q0(t), t)**2/2
	return e
def get_U():
	g = Symbol('g')
	l_com = Symbol('l_com')
	m = Symbol('m')
	t = Symbol('t')
	l_m = Symbol('l_m')
	m_1 = Symbol('m_1')
	q0 = Function('q0')
	e = g*l_com*m*sin(q0(t)) + g*l_m*m_1*sin(q0(t))
	return e
def get_Q():
	V = Symbol('V')
	c = Symbol('c')
	rho = Symbol('rho')
	s = Symbol('s')
	C_Dmax = Symbol('C_Dmax')
	alpha_1 = Symbol('alpha_1')
	t = Symbol('t')
	a0 = Symbol('a0')
	q0 = Function('q0')
	e = MutableDenseMatrix([[V**2*c*rho*s**2*(-C_Dmax*(1 - cos(2*alpha_1 - s*Derivative(q0(t), t)/V))*sin(alpha_1 - s*Derivative(q0(t), t)/(2*V))/2 + a0*(alpha_1 - s*Derivative(q0(t), t)/(2*V))*cos(alpha_1 - s*Derivative(q0(t), t)/(2*V)))/4]])
	return e
