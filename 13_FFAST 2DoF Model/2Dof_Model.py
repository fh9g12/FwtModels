from sympy import *
def get_M():
	m = Symbol('m')
	m_w = Symbol('m_w')
	l_com = Symbol('l_com')
	t = Symbol('t')
	I_xx = Symbol('I_xx')
	q1 = Function('q1')
	e = MutableDenseMatrix([[m + m_w, -l_com*m*cos(q1(t))], [-l_com*m*cos(q1(t)), I_xx + l_com**2*m]])
	return e
def get_f():
	g = Symbol('g')
	m = Symbol('m')
	m_w = Symbol('m_w')
	k_w = Symbol('k_w')
	t = Symbol('t')
	l_com = Symbol('l_com')
	k_fwt = Symbol('k_fwt')
	q0 = Function('q0')
	q1 = Function('q1')
	e = ImmutableDenseMatrix([[g*m + g*m_w + k_w*q0(t) + l_com*m*sin(q1(t))*Derivative(q1(t), t)**2], [-g*l_com*m*cos(q1(t)) + k_fwt*q1(t)]])
	return e
def get_T():
	I_xx = Symbol('I_xx')
	t = Symbol('t')
	l_com = Symbol('l_com')
	m = Symbol('m')
	m_w = Symbol('m_w')
	q1 = Function('q1')
	q0 = Function('q0')
	e = I_xx*Derivative(q1(t), t)**2/2 + l_com**2*m*Derivative(q1(t), t)**2/2 - l_com*m*cos(q1(t))*Derivative(q0(t), t)*Derivative(q1(t), t) + m*Derivative(q0(t), t)**2/2 + m_w*Derivative(q0(t), t)**2/2
	return e
def get_U():
	g = Symbol('g')
	m = Symbol('m')
	l_com = Symbol('l_com')
	t = Symbol('t')
	m_w = Symbol('m_w')
	k_fwt = Symbol('k_fwt')
	k_w = Symbol('k_w')
	q1 = Function('q1')
	q0 = Function('q0')
	e = g*m*(-l_com*sin(q1(t)) + q0(t)) + g*m_w*q0(t) + k_fwt*q1(t)**2/2 + k_w*q0(t)**2/2
	return e
def get_Q():
	V = Symbol('V')
	a_0 = Symbol('a_0')
	c = Symbol('c')
	rho = Symbol('rho')
	s = Symbol('s')
	alpha_1 = Symbol('alpha_1')
	t = Symbol('t')
	q1 = Function('q1')
	q0 = Function('q0')
	e = MutableDenseMatrix([[V**2*a_0*c*rho*s*(alpha_1 - (-s*Derivative(q1(t), t)/2 + cos(q1(t))*Derivative(q0(t), t))/V)*cos(q1(t))/2], [-V**2*a_0*c*rho*s**2*(alpha_1 - (-s*Derivative(q1(t), t)/2 + cos(q1(t))*Derivative(q0(t), t))/V)/4]])
	return e
