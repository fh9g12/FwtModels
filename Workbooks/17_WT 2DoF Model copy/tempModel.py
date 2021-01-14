from sympy import *
def get_M():
	m = Symbol('m')
	m_w = Symbol('m_w')
	l_com = Symbol('l_com')
	t = Symbol('t')
	q1 = Function('q1')
	e = MutableDenseMatrix([[m + m_w, -l_com*m*cos(q1(t))], [-l_com*m*cos(q1(t)), l_com**2*m]])
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
	m = Symbol('m')
	l_com = Symbol('l_com')
	t = Symbol('t')
	m_w = Symbol('m_w')
	q1 = Function('q1')
	q0 = Function('q0')
	e = m*(l_com**2*Derivative(q1(t), t)**2/2 - l_com*cos(q1(t))*Derivative(q0(t), t)*Derivative(q1(t), t) + Derivative(q0(t), t)**2/2) + m_w*Derivative(q0(t), t)**2/2
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
	c_root = Symbol('c_root')
	rho = Symbol('rho')
	s = Symbol('s')
	C_Dmax = Symbol('C_Dmax')
	Lambda = Symbol('Lambda')
	t = Symbol('t')
	w_g = Symbol('w_g')
	a0 = Symbol('a0')
	mu = Symbol('mu')
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
	e = MutableDenseMatrix([[V**2*c_root*rho*s*(-C_Dmax*(1 - cos(2*atan(sin(Lambda)*sin(q1(t))) - 2*w_g/V + 2*cos(q1(t))*Derivative(q0(t), t)/V))*sin(atan(sin(Lambda)*sin(q1(t))) + cos(q1(t))*Derivative(q0(t), t)/V)/2 + a0*(Float('-0.23999999999999999', precision=53) + log((exp(mu*(-atan(sin(Lambda)*sin(q1(t))) + Float('0.23999999999999999', precision=53) + w_g/V - cos(q1(t))*Derivative(q0(t), t)/V)) + 1)/(exp(mu*(-atan(sin(Lambda)*sin(q1(t))) + Float('-0.23999999999999999', precision=53) + w_g/V - cos(q1(t))*Derivative(q0(t), t)/V)) + 1))/mu)*cos(atan(sin(Lambda)*sin(q1(t))) + cos(q1(t))*Derivative(q0(t), t)/V))*cos(q1(t))/20 + V**2*c_root*rho*s*(-C_Dmax*(1 - cos(2*atan(sin(Lambda)*sin(q1(t))) - 2*w_g/V + 2*cos(q1(t))*Derivative(q0(t), t)/V))*sin(atan(sin(Lambda)*sin(q1(t))) + cos(q1(t))*Derivative(q0(t), t)/V)/2 + a1*(Float('-0.23999999999999999', precision=53) + log((exp(mu*(-atan(sin(Lambda)*sin(q1(t))) + Float('0.23999999999999999', precision=53) + w_g/V - cos(q1(t))*Derivative(q0(t), t)/V)) + 1)/(exp(mu*(-atan(sin(Lambda)*sin(q1(t))) + Float('-0.23999999999999999', precision=53) + w_g/V - cos(q1(t))*Derivative(q0(t), t)/V)) + 1))/mu)*cos(atan(sin(Lambda)*sin(q1(t))) + cos(q1(t))*Derivative(q0(t), t)/V))*cos(q1(t))/20 + V**2*c_root*rho*s*(-C_Dmax*(1 - cos(2*atan(sin(Lambda)*sin(q1(t))) - 2*w_g/V + 2*cos(q1(t))*Derivative(q0(t), t)/V))*sin(atan(sin(Lambda)*sin(q1(t))) + cos(q1(t))*Derivative(q0(t), t)/V)/2 + a2*(Float('-0.23999999999999999', precision=53) + log((exp(mu*(-atan(sin(Lambda)*sin(q1(t))) + Float('0.23999999999999999', precision=53) + w_g/V - cos(q1(t))*Derivative(q0(t), t)/V)) + 1)/(exp(mu*(-atan(sin(Lambda)*sin(q1(t))) + Float('-0.23999999999999999', precision=53) + w_g/V - cos(q1(t))*Derivative(q0(t), t)/V)) + 1))/mu)*cos(atan(sin(Lambda)*sin(q1(t))) + cos(q1(t))*Derivative(q0(t), t)/V))*cos(q1(t))/20 + V**2*c_root*rho*s*(-C_Dmax*(1 - cos(2*atan(sin(Lambda)*sin(q1(t))) - 2*w_g/V + 2*cos(q1(t))*Derivative(q0(t), t)/V))*sin(atan(sin(Lambda)*sin(q1(t))) + cos(q1(t))*Derivative(q0(t), t)/V)/2 + a3*(Float('-0.23999999999999999', precision=53) + log((exp(mu*(-atan(sin(Lambda)*sin(q1(t))) + Float('0.23999999999999999', precision=53) + w_g/V - cos(q1(t))*Derivative(q0(t), t)/V)) + 1)/(exp(mu*(-atan(sin(Lambda)*sin(q1(t))) + Float('-0.23999999999999999', precision=53) + w_g/V - cos(q1(t))*Derivative(q0(t), t)/V)) + 1))/mu)*cos(atan(sin(Lambda)*sin(q1(t))) + cos(q1(t))*Derivative(q0(t), t)/V))*cos(q1(t))/20 + V**2*c_root*rho*s*(-C_Dmax*(1 - cos(2*atan(sin(Lambda)*sin(q1(t))) - 2*w_g/V + 2*cos(q1(t))*Derivative(q0(t), t)/V))*sin(atan(sin(Lambda)*sin(q1(t))) + cos(q1(t))*Derivative(q0(t), t)/V)/2 + a4*(Float('-0.23999999999999999', precision=53) + log((exp(mu*(-atan(sin(Lambda)*sin(q1(t))) + Float('0.23999999999999999', precision=53) + w_g/V - cos(q1(t))*Derivative(q0(t), t)/V)) + 1)/(exp(mu*(-atan(sin(Lambda)*sin(q1(t))) + Float('-0.23999999999999999', precision=53) + w_g/V - cos(q1(t))*Derivative(q0(t), t)/V)) + 1))/mu)*cos(atan(sin(Lambda)*sin(q1(t))) + cos(q1(t))*Derivative(q0(t), t)/V))*cos(q1(t))/20 + V**2*c_root*rho*s*(-C_Dmax*(1 - cos(2*atan(sin(Lambda)*sin(q1(t))) - 2*w_g/V + 2*cos(q1(t))*Derivative(q0(t), t)/V))*sin(atan(sin(Lambda)*sin(q1(t))) + cos(q1(t))*Derivative(q0(t), t)/V)/2 + a5*(Float('-0.23999999999999999', precision=53) + log((exp(mu*(-atan(sin(Lambda)*sin(q1(t))) + Float('0.23999999999999999', precision=53) + w_g/V - cos(q1(t))*Derivative(q0(t), t)/V)) + 1)/(exp(mu*(-atan(sin(Lambda)*sin(q1(t))) + Float('-0.23999999999999999', precision=53) + w_g/V - cos(q1(t))*Derivative(q0(t), t)/V)) + 1))/mu)*cos(atan(sin(Lambda)*sin(q1(t))) + cos(q1(t))*Derivative(q0(t), t)/V))*cos(q1(t))/20 + V**2*c_root*rho*s*(-C_Dmax*(1 - cos(2*atan(sin(Lambda)*sin(q1(t))) - 2*w_g/V + 2*cos(q1(t))*Derivative(q0(t), t)/V))*sin(atan(sin(Lambda)*sin(q1(t))) + cos(q1(t))*Derivative(q0(t), t)/V)/2 + a6*(Float('-0.23999999999999999', precision=53) + log((exp(mu*(-atan(sin(Lambda)*sin(q1(t))) + Float('0.23999999999999999', precision=53) + w_g/V - cos(q1(t))*Derivative(q0(t), t)/V)) + 1)/(exp(mu*(-atan(sin(Lambda)*sin(q1(t))) + Float('-0.23999999999999999', precision=53) + w_g/V - cos(q1(t))*Derivative(q0(t), t)/V)) + 1))/mu)*cos(atan(sin(Lambda)*sin(q1(t))) + cos(q1(t))*Derivative(q0(t), t)/V))*cos(q1(t))/20 + V**2*c_root*rho*s*(-C_Dmax*(1 - cos(2*atan(sin(Lambda)*sin(q1(t))) - 2*w_g/V + 2*cos(q1(t))*Derivative(q0(t), t)/V))*sin(atan(sin(Lambda)*sin(q1(t))) + cos(q1(t))*Derivative(q0(t), t)/V)/2 + a7*(Float('-0.23999999999999999', precision=53) + log((exp(mu*(-atan(sin(Lambda)*sin(q1(t))) + Float('0.23999999999999999', precision=53) + w_g/V - cos(q1(t))*Derivative(q0(t), t)/V)) + 1)/(exp(mu*(-atan(sin(Lambda)*sin(q1(t))) + Float('-0.23999999999999999', precision=53) + w_g/V - cos(q1(t))*Derivative(q0(t), t)/V)) + 1))/mu)*cos(atan(sin(Lambda)*sin(q1(t))) + cos(q1(t))*Derivative(q0(t), t)/V))*cos(q1(t))/20 + V**2*c_root*rho*s*(-C_Dmax*(1 - cos(2*atan(sin(Lambda)*sin(q1(t))) - 2*w_g/V + 2*cos(q1(t))*Derivative(q0(t), t)/V))*sin(atan(sin(Lambda)*sin(q1(t))) + cos(q1(t))*Derivative(q0(t), t)/V)/2 + a8*(Float('-0.23999999999999999', precision=53) + log((exp(mu*(-atan(sin(Lambda)*sin(q1(t))) + Float('0.23999999999999999', precision=53) + w_g/V - cos(q1(t))*Derivative(q0(t), t)/V)) + 1)/(exp(mu*(-atan(sin(Lambda)*sin(q1(t))) + Float('-0.23999999999999999', precision=53) + w_g/V - cos(q1(t))*Derivative(q0(t), t)/V)) + 1))/mu)*cos(atan(sin(Lambda)*sin(q1(t))) + cos(q1(t))*Derivative(q0(t), t)/V))*cos(q1(t))/20 + V**2*c_root*rho*s*(-C_Dmax*(1 - cos(2*atan(sin(Lambda)*sin(q1(t))) - 2*w_g/V + 2*cos(q1(t))*Derivative(q0(t), t)/V))*sin(atan(sin(Lambda)*sin(q1(t))) + cos(q1(t))*Derivative(q0(t), t)/V)/2 + a9*(Float('-0.23999999999999999', precision=53) + log((exp(mu*(-atan(sin(Lambda)*sin(q1(t))) + Float('0.23999999999999999', precision=53) + w_g/V - cos(q1(t))*Derivative(q0(t), t)/V)) + 1)/(exp(mu*(-atan(sin(Lambda)*sin(q1(t))) + Float('-0.23999999999999999', precision=53) + w_g/V - cos(q1(t))*Derivative(q0(t), t)/V)) + 1))/mu)*cos(atan(sin(Lambda)*sin(q1(t))) + cos(q1(t))*Derivative(q0(t), t)/V))*cos(q1(t))/20], [0]])
	return e
