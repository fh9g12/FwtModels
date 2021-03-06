from sympy import *
def get_M():
	c = Symbol('c')
	delta_m = Symbol('delta_m')
	rho_t = Symbol('rho_t')
	s_t = Symbol('s_t')
	sigma = Symbol('sigma')
	t = Symbol('t')
	q1 = Function('q1')
	e = MutableDenseMatrix([[c*delta_m*rho_t*s_t*sigma, -c*delta_m*rho_t*s_t**2*sigma**2*cos(q1(t))/2], [-c*delta_m*rho_t*s_t**2*sigma**2*cos(q1(t))/2, c*delta_m*rho_t*s_t**3*sigma**3/3]])
	return e
def get_f():
	EI = Symbol('EI')
	t = Symbol('t')
	c = Symbol('c')
	delta_m = Symbol('delta_m')
	g = Symbol('g')
	rho_t = Symbol('rho_t')
	s_t = Symbol('s_t')
	sigma = Symbol('sigma')
	q0 = Function('q0')
	q1 = Function('q1')
	e = ImmutableDenseMatrix([[EI*q0(t) + c*delta_m*g*rho_t*s_t*sigma + c*delta_m*rho_t*s_t**2*sigma**2*sin(q1(t))*Derivative(q1(t), t)**2/2], [-c*delta_m*g*rho_t*s_t**2*sigma**2*cos(q1(t))/2]])
	return e
def get_T():
	c = Symbol('c')
	delta_m = Symbol('delta_m')
	rho_t = Symbol('rho_t')
	s_t = Symbol('s_t')
	sigma = Symbol('sigma')
	t = Symbol('t')
	q1 = Function('q1')
	q0 = Function('q0')
	e = c*delta_m*rho_t*s_t*sigma*(s_t**2*sigma**2*Derivative(q1(t), t)**2 - 3*s_t*sigma*cos(q1(t))*Derivative(q0(t), t)*Derivative(q1(t), t) + 3*Derivative(q0(t), t)**2)/6
	return e
def get_U():
	EI = Symbol('EI')
	t = Symbol('t')
	c = Symbol('c')
	delta_m = Symbol('delta_m')
	g = Symbol('g')
	rho_t = Symbol('rho_t')
	s_t = Symbol('s_t')
	sigma = Symbol('sigma')
	q0 = Function('q0')
	q1 = Function('q1')
	e = EI*q0(t)**2/2 + c*delta_m*g*rho_t*s_t*sigma*(-s_t*sigma*sin(q1(t))/2 + q0(t))
	return e
def get_Q():
	V = Symbol('V')
	a_0 = Symbol('a_0')
	c = Symbol('c')
	rho = Symbol('rho')
	s_t = Symbol('s_t')
	sigma = Symbol('sigma')
	Lambda = Symbol('Lambda')
	t = Symbol('t')
	alpha_r = Symbol('alpha_r')
	q1 = Function('q1')
	q0 = Function('q0')
	e = ImmutableDenseMatrix([[V**2*a_0*c*rho*s_t*sigma*(atan((sin(Lambda)*sin(q1(t))*cos(alpha_r) + sin(alpha_r)*cos(q1(t)))/((sin(Lambda)**2*cos(q1(t)) + cos(Lambda)**2)*cos(alpha_r) - sin(Lambda)*sin(alpha_r)*sin(q1(t)))) - (-9*s_t*sigma*Derivative(q1(t), t)/10 + cos(q1(t))*Derivative(q0(t), t))/V)*cos(atan((sin(Lambda)*sin(q1(t))*cos(alpha_r) + sin(alpha_r)*cos(q1(t)))/((sin(Lambda)**2*cos(q1(t)) + cos(Lambda)**2)*cos(alpha_r) - sin(Lambda)*sin(alpha_r)*sin(q1(t)))) - (-9*s_t*sigma*Derivative(q1(t), t)/10 + cos(q1(t))*Derivative(q0(t), t))/V)*cos(q1(t))/100 + 3*V**2*a_0*c*rho*s_t*sigma*(atan((sin(Lambda)*sin(q1(t))*cos(alpha_r) + sin(alpha_r)*cos(q1(t)))/((sin(Lambda)**2*cos(q1(t)) + cos(Lambda)**2)*cos(alpha_r) - sin(Lambda)*sin(alpha_r)*sin(q1(t)))) - (-7*s_t*sigma*Derivative(q1(t), t)/10 + cos(q1(t))*Derivative(q0(t), t))/V)*cos(atan((sin(Lambda)*sin(q1(t))*cos(alpha_r) + sin(alpha_r)*cos(q1(t)))/((sin(Lambda)**2*cos(q1(t)) + cos(Lambda)**2)*cos(alpha_r) - sin(Lambda)*sin(alpha_r)*sin(q1(t)))) - (-7*s_t*sigma*Derivative(q1(t), t)/10 + cos(q1(t))*Derivative(q0(t), t))/V)*cos(q1(t))/100 + V**2*a_0*c*rho*s_t*sigma*(atan((sin(Lambda)*sin(q1(t))*cos(alpha_r) + sin(alpha_r)*cos(q1(t)))/((sin(Lambda)**2*cos(q1(t)) + cos(Lambda)**2)*cos(alpha_r) - sin(Lambda)*sin(alpha_r)*sin(q1(t)))) - (-s_t*sigma*Derivative(q1(t), t)/2 + cos(q1(t))*Derivative(q0(t), t))/V)*cos(atan((sin(Lambda)*sin(q1(t))*cos(alpha_r) + sin(alpha_r)*cos(q1(t)))/((sin(Lambda)**2*cos(q1(t)) + cos(Lambda)**2)*cos(alpha_r) - sin(Lambda)*sin(alpha_r)*sin(q1(t)))) - (-s_t*sigma*Derivative(q1(t), t)/2 + cos(q1(t))*Derivative(q0(t), t))/V)*cos(q1(t))/20 + 7*V**2*a_0*c*rho*s_t*sigma*(atan((sin(Lambda)*sin(q1(t))*cos(alpha_r) + sin(alpha_r)*cos(q1(t)))/((sin(Lambda)**2*cos(q1(t)) + cos(Lambda)**2)*cos(alpha_r) - sin(Lambda)*sin(alpha_r)*sin(q1(t)))) - (-3*s_t*sigma*Derivative(q1(t), t)/10 + cos(q1(t))*Derivative(q0(t), t))/V)*cos(atan((sin(Lambda)*sin(q1(t))*cos(alpha_r) + sin(alpha_r)*cos(q1(t)))/((sin(Lambda)**2*cos(q1(t)) + cos(Lambda)**2)*cos(alpha_r) - sin(Lambda)*sin(alpha_r)*sin(q1(t)))) - (-3*s_t*sigma*Derivative(q1(t), t)/10 + cos(q1(t))*Derivative(q0(t), t))/V)*cos(q1(t))/100 + 9*V**2*a_0*c*rho*s_t*sigma*(atan((sin(Lambda)*sin(q1(t))*cos(alpha_r) + sin(alpha_r)*cos(q1(t)))/((sin(Lambda)**2*cos(q1(t)) + cos(Lambda)**2)*cos(alpha_r) - sin(Lambda)*sin(alpha_r)*sin(q1(t)))) - (-s_t*sigma*Derivative(q1(t), t)/10 + cos(q1(t))*Derivative(q0(t), t))/V)*cos(atan((sin(Lambda)*sin(q1(t))*cos(alpha_r) + sin(alpha_r)*cos(q1(t)))/((sin(Lambda)**2*cos(q1(t)) + cos(Lambda)**2)*cos(alpha_r) - sin(Lambda)*sin(alpha_r)*sin(q1(t)))) - (-s_t*sigma*Derivative(q1(t), t)/10 + cos(q1(t))*Derivative(q0(t), t))/V)*cos(q1(t))/100], [-9*V**2*a_0*c*rho*s_t**2*sigma**2*(atan((sin(Lambda)*sin(q1(t))*cos(alpha_r) + sin(alpha_r)*cos(q1(t)))/((sin(Lambda)**2*cos(q1(t)) + cos(Lambda)**2)*cos(alpha_r) - sin(Lambda)*sin(alpha_r)*sin(q1(t)))) - (-9*s_t*sigma*Derivative(q1(t), t)/10 + cos(q1(t))*Derivative(q0(t), t))/V)*cos(atan((sin(Lambda)*sin(q1(t))*cos(alpha_r) + sin(alpha_r)*cos(q1(t)))/((sin(Lambda)**2*cos(q1(t)) + cos(Lambda)**2)*cos(alpha_r) - sin(Lambda)*sin(alpha_r)*sin(q1(t)))) - (-9*s_t*sigma*Derivative(q1(t), t)/10 + cos(q1(t))*Derivative(q0(t), t))/V)/1000 - 21*V**2*a_0*c*rho*s_t**2*sigma**2*(atan((sin(Lambda)*sin(q1(t))*cos(alpha_r) + sin(alpha_r)*cos(q1(t)))/((sin(Lambda)**2*cos(q1(t)) + cos(Lambda)**2)*cos(alpha_r) - sin(Lambda)*sin(alpha_r)*sin(q1(t)))) - (-7*s_t*sigma*Derivative(q1(t), t)/10 + cos(q1(t))*Derivative(q0(t), t))/V)*cos(atan((sin(Lambda)*sin(q1(t))*cos(alpha_r) + sin(alpha_r)*cos(q1(t)))/((sin(Lambda)**2*cos(q1(t)) + cos(Lambda)**2)*cos(alpha_r) - sin(Lambda)*sin(alpha_r)*sin(q1(t)))) - (-7*s_t*sigma*Derivative(q1(t), t)/10 + cos(q1(t))*Derivative(q0(t), t))/V)/1000 - V**2*a_0*c*rho*s_t**2*sigma**2*(atan((sin(Lambda)*sin(q1(t))*cos(alpha_r) + sin(alpha_r)*cos(q1(t)))/((sin(Lambda)**2*cos(q1(t)) + cos(Lambda)**2)*cos(alpha_r) - sin(Lambda)*sin(alpha_r)*sin(q1(t)))) - (-s_t*sigma*Derivative(q1(t), t)/2 + cos(q1(t))*Derivative(q0(t), t))/V)*cos(atan((sin(Lambda)*sin(q1(t))*cos(alpha_r) + sin(alpha_r)*cos(q1(t)))/((sin(Lambda)**2*cos(q1(t)) + cos(Lambda)**2)*cos(alpha_r) - sin(Lambda)*sin(alpha_r)*sin(q1(t)))) - (-s_t*sigma*Derivative(q1(t), t)/2 + cos(q1(t))*Derivative(q0(t), t))/V)/40 - 21*V**2*a_0*c*rho*s_t**2*sigma**2*(atan((sin(Lambda)*sin(q1(t))*cos(alpha_r) + sin(alpha_r)*cos(q1(t)))/((sin(Lambda)**2*cos(q1(t)) + cos(Lambda)**2)*cos(alpha_r) - sin(Lambda)*sin(alpha_r)*sin(q1(t)))) - (-3*s_t*sigma*Derivative(q1(t), t)/10 + cos(q1(t))*Derivative(q0(t), t))/V)*cos(atan((sin(Lambda)*sin(q1(t))*cos(alpha_r) + sin(alpha_r)*cos(q1(t)))/((sin(Lambda)**2*cos(q1(t)) + cos(Lambda)**2)*cos(alpha_r) - sin(Lambda)*sin(alpha_r)*sin(q1(t)))) - (-3*s_t*sigma*Derivative(q1(t), t)/10 + cos(q1(t))*Derivative(q0(t), t))/V)/1000 - 9*V**2*a_0*c*rho*s_t**2*sigma**2*(atan((sin(Lambda)*sin(q1(t))*cos(alpha_r) + sin(alpha_r)*cos(q1(t)))/((sin(Lambda)**2*cos(q1(t)) + cos(Lambda)**2)*cos(alpha_r) - sin(Lambda)*sin(alpha_r)*sin(q1(t)))) - (-s_t*sigma*Derivative(q1(t), t)/10 + cos(q1(t))*Derivative(q0(t), t))/V)*cos(atan((sin(Lambda)*sin(q1(t))*cos(alpha_r) + sin(alpha_r)*cos(q1(t)))/((sin(Lambda)**2*cos(q1(t)) + cos(Lambda)**2)*cos(alpha_r) - sin(Lambda)*sin(alpha_r)*sin(q1(t)))) - (-s_t*sigma*Derivative(q1(t), t)/10 + cos(q1(t))*Derivative(q0(t), t))/V)/1000]])
	return e
