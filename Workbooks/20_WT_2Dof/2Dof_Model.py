from sympy import *
def get_M():
	m = Symbol('m')
	m_1 = Symbol('m_1')
	m_w = Symbol('m_w')
	l_com = Symbol('l_com')
	t = Symbol('t')
	l_m = Symbol('l_m')
	q1 = Function('q1')
	e = MutableDenseMatrix([[m + m_1 + m_w, l_com*m*cos(q1(t)) + l_m*m_1*cos(q1(t))], [l_com*m*cos(q1(t)) + l_m*m_1*cos(q1(t)), l_com**2*m + l_m**2*m_1]])
	return e
def get_f():
	d_w = Symbol('d_w')
	t = Symbol('t')
	g = Symbol('g')
	m = Symbol('m')
	m_1 = Symbol('m_1')
	m_w = Symbol('m_w')
	k_w = Symbol('k_w')
	l_com = Symbol('l_com')
	l_m = Symbol('l_m')
	q0 = Function('q0')
	q1 = Function('q1')
	e = ImmutableDenseMatrix([[d_w*Derivative(q0(t), t) + g*m + g*m_1 + g*m_w + k_w*q0(t) - l_com*m*sin(q1(t))*Derivative(q1(t), t)**2 - l_m*m_1*sin(q1(t))*Derivative(q1(t), t)**2], [g*l_com*m*cos(q1(t)) + g*l_m*m_1*cos(q1(t))]])
	return e
def get_T():
	m = Symbol('m')
	l_com = Symbol('l_com')
	t = Symbol('t')
	m_1 = Symbol('m_1')
	l_m = Symbol('l_m')
	m_w = Symbol('m_w')
	q1 = Function('q1')
	q0 = Function('q0')
	e = m*(l_com**2*Derivative(q1(t), t)**2/2 + l_com*cos(q1(t))*Derivative(q0(t), t)*Derivative(q1(t), t) + Derivative(q0(t), t)**2/2) + m_1*(l_m**2*Derivative(q1(t), t)**2/2 + l_m*cos(q1(t))*Derivative(q0(t), t)*Derivative(q1(t), t) + Derivative(q0(t), t)**2/2) + m_w*Derivative(q0(t), t)**2/2
	return e
def get_U():
	g = Symbol('g')
	m = Symbol('m')
	l_com = Symbol('l_com')
	t = Symbol('t')
	m_1 = Symbol('m_1')
	l_m = Symbol('l_m')
	m_w = Symbol('m_w')
	k_w = Symbol('k_w')
	q1 = Function('q1')
	q0 = Function('q0')
	e = g*m*(l_com*sin(q1(t)) + q0(t)) + g*m_1*(l_m*sin(q1(t)) + q0(t)) + g*m_w*q0(t) + k_w*q0(t)**2/2
	return e
def get_Q():
	V = Symbol('V')
	c = Symbol('c')
	rho = Symbol('rho')
	s = Symbol('s')
	a0 = Symbol('a0')
	alpha_1 = Symbol('alpha_1')
	t = Symbol('t')
	d_a = Symbol('d_a')
	q1 = Function('q1')
	q0 = Function('q0')
	e = MutableDenseMatrix([[V**2*c*rho*s*a0*(alpha_1 - (s*Derivative(q1(t), t)/2 + cos(q1(t))*Derivative(q0(t), t))/V)*cos(q1(t))/2 - V*c*d_a*rho*Derivative(q0(t), t)/2], [V**2*c*rho*s**2*a0*(alpha_1 - (s*Derivative(q1(t), t)/2 + cos(q1(t))*Derivative(q0(t), t))/V)/4]])
	return e
