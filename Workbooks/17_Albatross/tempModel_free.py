from sympy import *
def get_M():
	l_com = Symbol('l_com')
	m = Symbol('m')
	e = MutableDenseMatrix([[l_com**2*m]])
	return e
def get_f():
	g = Symbol('g')
	l_com = Symbol('l_com')
	m = Symbol('m')
	alpha_r = Symbol('alpha_r')
	t = Symbol('t')
	q0 = Function('q0')
	e = ImmutableDenseMatrix([[g*l_com*m*cos(alpha_r)*cos(q0(t))]])
	return e
def get_T():
	l_com = Symbol('l_com')
	m = Symbol('m')
	t = Symbol('t')
	q0 = Function('q0')
	e = l_com**2*m*Derivative(q0(t), t)**2/2
	return e
def get_U():
	g = Symbol('g')
	l_com = Symbol('l_com')
	m = Symbol('m')
	t = Symbol('t')
	alpha_r = Symbol('alpha_r')
	q0 = Function('q0')
	e = g*l_com*m*sin(q0(t))*cos(alpha_r)
	return e
def get_Q():
	V = Symbol('V')
	c = Symbol('c')
	rho = Symbol('rho')
	s = Symbol('s')
	a0 = Symbol('a0')
	alpha_c = Symbol('alpha_c')
	alpha_r = Symbol('alpha_r')
	t = Symbol('t')
	Lambda = Symbol('Lambda')
	beta = Symbol('beta')
	a1 = Symbol('a1')
	a2 = Symbol('a2')
	a3 = Symbol('a3')
	a4 = Symbol('a4')
	a5 = Symbol('a5')
	a6 = Symbol('a6')
	a7 = Symbol('a7')
	a8 = Symbol('a8')
	a9 = Symbol('a9')
	q0 = Function('q0')
	e = MutableDenseMatrix([[V**2*c*rho*s**2*a0*(alpha_c - atan((sin(alpha_r)*cos(q0(t)) + sin(Lambda + beta)*sin(q0(t))*cos(alpha_r))/(-sin(Lambda)*sin(alpha_r)*sin(q0(t)) + sin(Lambda)*sin(Lambda + beta)*cos(alpha_r)*cos(q0(t)) - sin(Lambda)*sin(Lambda + beta)*cos(alpha_r) + cos(alpha_r)*cos(beta))) - s*Derivative(q0(t), t)/(20*V*cos(beta)))*cos(beta)**2/400 + 3*V**2*c*rho*s**2*a1*(alpha_c - atan((sin(alpha_r)*cos(q0(t)) + sin(Lambda + beta)*sin(q0(t))*cos(alpha_r))/(-sin(Lambda)*sin(alpha_r)*sin(q0(t)) + sin(Lambda)*sin(Lambda + beta)*cos(alpha_r)*cos(q0(t)) - sin(Lambda)*sin(Lambda + beta)*cos(alpha_r) + cos(alpha_r)*cos(beta))) - 3*s*Derivative(q0(t), t)/(20*V*cos(beta)))*cos(beta)**2/400 + V**2*c*rho*s**2*a2*(alpha_c - atan((sin(alpha_r)*cos(q0(t)) + sin(Lambda + beta)*sin(q0(t))*cos(alpha_r))/(-sin(Lambda)*sin(alpha_r)*sin(q0(t)) + sin(Lambda)*sin(Lambda + beta)*cos(alpha_r)*cos(q0(t)) - sin(Lambda)*sin(Lambda + beta)*cos(alpha_r) + cos(alpha_r)*cos(beta))) - s*Derivative(q0(t), t)/(4*V*cos(beta)))*cos(beta)**2/80 + 7*V**2*c*rho*s**2*a3*(alpha_c - atan((sin(alpha_r)*cos(q0(t)) + sin(Lambda + beta)*sin(q0(t))*cos(alpha_r))/(-sin(Lambda)*sin(alpha_r)*sin(q0(t)) + sin(Lambda)*sin(Lambda + beta)*cos(alpha_r)*cos(q0(t)) - sin(Lambda)*sin(Lambda + beta)*cos(alpha_r) + cos(alpha_r)*cos(beta))) - 7*s*Derivative(q0(t), t)/(20*V*cos(beta)))*cos(beta)**2/400 + 9*V**2*c*rho*s**2*a4*(alpha_c - atan((sin(alpha_r)*cos(q0(t)) + sin(Lambda + beta)*sin(q0(t))*cos(alpha_r))/(-sin(Lambda)*sin(alpha_r)*sin(q0(t)) + sin(Lambda)*sin(Lambda + beta)*cos(alpha_r)*cos(q0(t)) - sin(Lambda)*sin(Lambda + beta)*cos(alpha_r) + cos(alpha_r)*cos(beta))) - 9*s*Derivative(q0(t), t)/(20*V*cos(beta)))*cos(beta)**2/400 + 11*V**2*c*rho*s**2*a5*(alpha_c - atan((sin(alpha_r)*cos(q0(t)) + sin(Lambda + beta)*sin(q0(t))*cos(alpha_r))/(-sin(Lambda)*sin(alpha_r)*sin(q0(t)) + sin(Lambda)*sin(Lambda + beta)*cos(alpha_r)*cos(q0(t)) - sin(Lambda)*sin(Lambda + beta)*cos(alpha_r) + cos(alpha_r)*cos(beta))) - 11*s*Derivative(q0(t), t)/(20*V*cos(beta)))*cos(beta)**2/400 + 13*V**2*c*rho*s**2*a6*(alpha_c - atan((sin(alpha_r)*cos(q0(t)) + sin(Lambda + beta)*sin(q0(t))*cos(alpha_r))/(-sin(Lambda)*sin(alpha_r)*sin(q0(t)) + sin(Lambda)*sin(Lambda + beta)*cos(alpha_r)*cos(q0(t)) - sin(Lambda)*sin(Lambda + beta)*cos(alpha_r) + cos(alpha_r)*cos(beta))) - 13*s*Derivative(q0(t), t)/(20*V*cos(beta)))*cos(beta)**2/400 + 3*V**2*c*rho*s**2*a7*(alpha_c - atan((sin(alpha_r)*cos(q0(t)) + sin(Lambda + beta)*sin(q0(t))*cos(alpha_r))/(-sin(Lambda)*sin(alpha_r)*sin(q0(t)) + sin(Lambda)*sin(Lambda + beta)*cos(alpha_r)*cos(q0(t)) - sin(Lambda)*sin(Lambda + beta)*cos(alpha_r) + cos(alpha_r)*cos(beta))) - 3*s*Derivative(q0(t), t)/(4*V*cos(beta)))*cos(beta)**2/80 + 17*V**2*c*rho*s**2*a8*(alpha_c - atan((sin(alpha_r)*cos(q0(t)) + sin(Lambda + beta)*sin(q0(t))*cos(alpha_r))/(-sin(Lambda)*sin(alpha_r)*sin(q0(t)) + sin(Lambda)*sin(Lambda + beta)*cos(alpha_r)*cos(q0(t)) - sin(Lambda)*sin(Lambda + beta)*cos(alpha_r) + cos(alpha_r)*cos(beta))) - 17*s*Derivative(q0(t), t)/(20*V*cos(beta)))*cos(beta)**2/400 + 19*V**2*c*rho*s**2*a9*(alpha_c - atan((sin(alpha_r)*cos(q0(t)) + sin(Lambda + beta)*sin(q0(t))*cos(alpha_r))/(-sin(Lambda)*sin(alpha_r)*sin(q0(t)) + sin(Lambda)*sin(Lambda + beta)*cos(alpha_r)*cos(q0(t)) - sin(Lambda)*sin(Lambda + beta)*cos(alpha_r) + cos(alpha_r)*cos(beta))) - 19*s*Derivative(q0(t), t)/(20*V*cos(beta)))*cos(beta)**2/400]])
	return e
