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
	q0 = Function('q0')
	q1 = Function('q1')
	e = ImmutableDenseMatrix([[g*m + g*m_w + k_w*q0(t) + l_com*m*sin(q1(t))*Derivative(q1(t), t)**2], [-g*l_com*m*cos(q1(t))]])
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
	k_w = Symbol('k_w')
	q1 = Function('q1')
	q0 = Function('q0')
	e = g*m*(-l_com*sin(q1(t)) + q0(t)) + g*m_w*q0(t) + k_w*q0(t)**2/2
	return e
def get_Q():
	V = Symbol('V')
	a0 = Symbol('a0')
	c = Symbol('c')
	rho = Symbol('rho')
	s = Symbol('s')
	alpha_r = Symbol('alpha_r')
	t = Symbol('t')
	Lambda = Symbol('Lambda')
	a1 = Symbol('a1')
	a2 = Symbol('a2')
	a3 = Symbol('a3')
	a4 = Symbol('a4')
	a5 = Symbol('a5')
	a6 = Symbol('a6')
	a7 = Symbol('a7')
	a8 = Symbol('a8')
	a9 = Symbol('a9')
	q1 = Function('q1')
	q0 = Function('q0')
	e = MutableDenseMatrix([[V**2*a0*c*rho*s*(alpha_r*cos(q1(t)) + atan(sin(Lambda)*sin(q1(t))) - (-s*Derivative(q1(t), t)/20 + cos(q1(t))*Derivative(q0(t), t))/V)*cos(q1(t))/20 + V**2*a1*c*rho*s*(alpha_r*cos(q1(t)) + atan(sin(Lambda)*sin(q1(t))) - (-3*s*Derivative(q1(t), t)/20 + cos(q1(t))*Derivative(q0(t), t))/V)*cos(q1(t))/20 + V**2*a2*c*rho*s*(alpha_r*cos(q1(t)) + atan(sin(Lambda)*sin(q1(t))) - (-s*Derivative(q1(t), t)/4 + cos(q1(t))*Derivative(q0(t), t))/V)*cos(q1(t))/20 + V**2*a3*c*rho*s*(alpha_r*cos(q1(t)) + atan(sin(Lambda)*sin(q1(t))) - (-7*s*Derivative(q1(t), t)/20 + cos(q1(t))*Derivative(q0(t), t))/V)*cos(q1(t))/20 + V**2*a4*c*rho*s*(alpha_r*cos(q1(t)) + atan(sin(Lambda)*sin(q1(t))) - (-9*s*Derivative(q1(t), t)/20 + cos(q1(t))*Derivative(q0(t), t))/V)*cos(q1(t))/20 + V**2*a5*c*rho*s*(alpha_r*cos(q1(t)) + atan(sin(Lambda)*sin(q1(t))) - (-11*s*Derivative(q1(t), t)/20 + cos(q1(t))*Derivative(q0(t), t))/V)*cos(q1(t))/20 + V**2*a6*c*rho*s*(alpha_r*cos(q1(t)) + atan(sin(Lambda)*sin(q1(t))) - (-13*s*Derivative(q1(t), t)/20 + cos(q1(t))*Derivative(q0(t), t))/V)*cos(q1(t))/20 + V**2*a7*c*rho*s*(alpha_r*cos(q1(t)) + atan(sin(Lambda)*sin(q1(t))) - (-3*s*Derivative(q1(t), t)/4 + cos(q1(t))*Derivative(q0(t), t))/V)*cos(q1(t))/20 + V**2*a8*c*rho*s*(alpha_r*cos(q1(t)) + atan(sin(Lambda)*sin(q1(t))) - (-17*s*Derivative(q1(t), t)/20 + cos(q1(t))*Derivative(q0(t), t))/V)*cos(q1(t))/20 + V**2*a9*c*rho*s*(alpha_r*cos(q1(t)) + atan(sin(Lambda)*sin(q1(t))) - (-19*s*Derivative(q1(t), t)/20 + cos(q1(t))*Derivative(q0(t), t))/V)*cos(q1(t))/20], [-V**2*a0*c*rho*s**2*(alpha_r*cos(q1(t)) + atan(sin(Lambda)*sin(q1(t))) - (-s*Derivative(q1(t), t)/20 + cos(q1(t))*Derivative(q0(t), t))/V)/400 - 3*V**2*a1*c*rho*s**2*(alpha_r*cos(q1(t)) + atan(sin(Lambda)*sin(q1(t))) - (-3*s*Derivative(q1(t), t)/20 + cos(q1(t))*Derivative(q0(t), t))/V)/400 - V**2*a2*c*rho*s**2*(alpha_r*cos(q1(t)) + atan(sin(Lambda)*sin(q1(t))) - (-s*Derivative(q1(t), t)/4 + cos(q1(t))*Derivative(q0(t), t))/V)/80 - 7*V**2*a3*c*rho*s**2*(alpha_r*cos(q1(t)) + atan(sin(Lambda)*sin(q1(t))) - (-7*s*Derivative(q1(t), t)/20 + cos(q1(t))*Derivative(q0(t), t))/V)/400 - 9*V**2*a4*c*rho*s**2*(alpha_r*cos(q1(t)) + atan(sin(Lambda)*sin(q1(t))) - (-9*s*Derivative(q1(t), t)/20 + cos(q1(t))*Derivative(q0(t), t))/V)/400 - 11*V**2*a5*c*rho*s**2*(alpha_r*cos(q1(t)) + atan(sin(Lambda)*sin(q1(t))) - (-11*s*Derivative(q1(t), t)/20 + cos(q1(t))*Derivative(q0(t), t))/V)/400 - 13*V**2*a6*c*rho*s**2*(alpha_r*cos(q1(t)) + atan(sin(Lambda)*sin(q1(t))) - (-13*s*Derivative(q1(t), t)/20 + cos(q1(t))*Derivative(q0(t), t))/V)/400 - 3*V**2*a7*c*rho*s**2*(alpha_r*cos(q1(t)) + atan(sin(Lambda)*sin(q1(t))) - (-3*s*Derivative(q1(t), t)/4 + cos(q1(t))*Derivative(q0(t), t))/V)/80 - 17*V**2*a8*c*rho*s**2*(alpha_r*cos(q1(t)) + atan(sin(Lambda)*sin(q1(t))) - (-17*s*Derivative(q1(t), t)/20 + cos(q1(t))*Derivative(q0(t), t))/V)/400 - 19*V**2*a9*c*rho*s**2*(alpha_r*cos(q1(t)) + atan(sin(Lambda)*sin(q1(t))) - (-19*s*Derivative(q1(t), t)/20 + cos(q1(t))*Derivative(q0(t), t))/V)/400]])
	return e
