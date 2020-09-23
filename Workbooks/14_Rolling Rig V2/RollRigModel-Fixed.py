from sympy import *
def get_M():
	m_f = Symbol('m_f')
	s = Symbol('s')
	sigma = Symbol('sigma')
	m_l = Symbol('m_l')
	m_w = Symbol('m_w')
	e = MutableDenseMatrix([[m_f*s**2*sigma**2/6 - m_f*s**2*sigma/2 + m_f*s**2/2 + m_l*s**2*sigma**2/2 - m_l*s**2*sigma + m_l*s**2/2 + m_w*s**2*sigma**2/12 - m_w*s**2*sigma/6 + m_w*s**2/12]])
	return e
def get_f():
	g = Symbol('g')
	m_f = Symbol('m_f')
	s = Symbol('s')
	sigma = Symbol('sigma')
	t = Symbol('t')
	q0 = Function('q0')
	e = ImmutableDenseMatrix([[-g*m_f*(-s*sigma*cos(q0(t))/4 - s*(Rational(1, 2) - sigma/2)*cos(q0(t))) - g*m_f*(s*sigma*cos(q0(t))/4 - s*(sigma/2 + Rational(-1, 2))*cos(q0(t)))]])
	return e
def get_T():
	m_f = Symbol('m_f')
	s = Symbol('s')
	sigma = Symbol('sigma')
	t = Symbol('t')
	m_l = Symbol('m_l')
	m_w = Symbol('m_w')
	q0 = Function('q0')
	e = m_f*s**2*sigma**2*Derivative(q0(t), t)**2/12 - m_f*s**2*sigma*Derivative(q0(t), t)**2/4 + m_f*s**2*Derivative(q0(t), t)**2/4 + m_l*s**2*sigma**2*Derivative(q0(t), t)**2/4 - m_l*s**2*sigma*Derivative(q0(t), t)**2/2 + m_l*s**2*Derivative(q0(t), t)**2/4 + m_w*s**2*sigma**2*Derivative(q0(t), t)**2/24 - m_w*s**2*sigma*Derivative(q0(t), t)**2/12 + m_w*s**2*Derivative(q0(t), t)**2/24
	return e
def get_U():
	g = Symbol('g')
	m_f = Symbol('m_f')
	s = Symbol('s')
	sigma = Symbol('sigma')
	t = Symbol('t')
	q0 = Function('q0')
	e = g*m_f*(-s*sigma*sin(q0(t))/4 - s*(1 - sigma)*sin(q0(t))/2) + g*m_f*(s*sigma*sin(q0(t))/4 + s*(1 - sigma)*sin(q0(t))/2)
	return e
def get_Q():
	V = Symbol('V')
	c = Symbol('c')
	rho = Symbol('rho')
	s = Symbol('s')
	sigma = Symbol('sigma')
	C_Dmax = Symbol('C_Dmax')
	alpha_1 = Symbol('alpha_1')
	w_g = Symbol('w_g')
	t = Symbol('t')
	a_1 = Symbol('a_1')
	alpha_2 = Symbol('alpha_2')
	a_0 = Symbol('a_0')
	q0 = Function('q0')
	e = ImmutableDenseMatrix([[V**2*c*rho*s*sigma*(-9*s*sigma/20 + s/2)*(C_Dmax*(1 - cos(2*alpha_1 + 2*w_g/V - 2*(-9*s*sigma/20 + s/2)*Derivative(q0(t), t)/V))*sin(alpha_1 - (-9*s*sigma/20 + s/2)*Derivative(q0(t), t)/V)/2 + a_1*(alpha_1 + w_g/V - (-9*s*sigma/20 + s/2)*Derivative(q0(t), t)/V)*cos(alpha_1 - (-9*s*sigma/20 + s/2)*Derivative(q0(t), t)/V))/20 + V**2*c*rho*s*sigma*(-7*s*sigma/20 + s/2)*(C_Dmax*(1 - cos(2*alpha_1 + 2*w_g/V - 2*(-7*s*sigma/20 + s/2)*Derivative(q0(t), t)/V))*sin(alpha_1 - (-7*s*sigma/20 + s/2)*Derivative(q0(t), t)/V)/2 + a_1*(alpha_1 + w_g/V - (-7*s*sigma/20 + s/2)*Derivative(q0(t), t)/V)*cos(alpha_1 - (-7*s*sigma/20 + s/2)*Derivative(q0(t), t)/V))/20 + V**2*c*rho*s*sigma*(-s*sigma/4 + s/2)*(C_Dmax*(1 - cos(2*alpha_1 + 2*w_g/V - 2*(-s*sigma/4 + s/2)*Derivative(q0(t), t)/V))*sin(alpha_1 - (-s*sigma/4 + s/2)*Derivative(q0(t), t)/V)/2 + a_1*(alpha_1 + w_g/V - (-s*sigma/4 + s/2)*Derivative(q0(t), t)/V)*cos(alpha_1 - (-s*sigma/4 + s/2)*Derivative(q0(t), t)/V))/20 + V**2*c*rho*s*sigma*(-3*s*sigma/20 + s/2)*(C_Dmax*(1 - cos(2*alpha_1 + 2*w_g/V - 2*(-3*s*sigma/20 + s/2)*Derivative(q0(t), t)/V))*sin(alpha_1 - (-3*s*sigma/20 + s/2)*Derivative(q0(t), t)/V)/2 + a_1*(alpha_1 + w_g/V - (-3*s*sigma/20 + s/2)*Derivative(q0(t), t)/V)*cos(alpha_1 - (-3*s*sigma/20 + s/2)*Derivative(q0(t), t)/V))/20 + V**2*c*rho*s*sigma*(-s*sigma/20 + s/2)*(C_Dmax*(1 - cos(2*alpha_1 + 2*w_g/V - 2*(-s*sigma/20 + s/2)*Derivative(q0(t), t)/V))*sin(alpha_1 - (-s*sigma/20 + s/2)*Derivative(q0(t), t)/V)/2 + a_1*(alpha_1 + w_g/V - (-s*sigma/20 + s/2)*Derivative(q0(t), t)/V)*cos(alpha_1 - (-s*sigma/20 + s/2)*Derivative(q0(t), t)/V))/20 + V**2*c*rho*s*sigma*(s*sigma/20 - s/2)*(C_Dmax*(1 - cos(2*alpha_2 + 2*w_g/V - 2*(s*sigma/20 - s/2)*Derivative(q0(t), t)/V))*sin(alpha_2 - (s*sigma/20 - s/2)*Derivative(q0(t), t)/V)/2 + a_1*(alpha_2 + w_g/V - (s*sigma/20 - s/2)*Derivative(q0(t), t)/V)*cos(alpha_2 - (s*sigma/20 - s/2)*Derivative(q0(t), t)/V))/20 + V**2*c*rho*s*sigma*(3*s*sigma/20 - s/2)*(C_Dmax*(1 - cos(2*alpha_2 + 2*w_g/V - 2*(3*s*sigma/20 - s/2)*Derivative(q0(t), t)/V))*sin(alpha_2 - (3*s*sigma/20 - s/2)*Derivative(q0(t), t)/V)/2 + a_1*(alpha_2 + w_g/V - (3*s*sigma/20 - s/2)*Derivative(q0(t), t)/V)*cos(alpha_2 - (3*s*sigma/20 - s/2)*Derivative(q0(t), t)/V))/20 + V**2*c*rho*s*sigma*(s*sigma/4 - s/2)*(C_Dmax*(1 - cos(2*alpha_2 + 2*w_g/V - 2*(s*sigma/4 - s/2)*Derivative(q0(t), t)/V))*sin(alpha_2 - (s*sigma/4 - s/2)*Derivative(q0(t), t)/V)/2 + a_1*(alpha_2 + w_g/V - (s*sigma/4 - s/2)*Derivative(q0(t), t)/V)*cos(alpha_2 - (s*sigma/4 - s/2)*Derivative(q0(t), t)/V))/20 + V**2*c*rho*s*sigma*(7*s*sigma/20 - s/2)*(C_Dmax*(1 - cos(2*alpha_2 + 2*w_g/V - 2*(7*s*sigma/20 - s/2)*Derivative(q0(t), t)/V))*sin(alpha_2 - (7*s*sigma/20 - s/2)*Derivative(q0(t), t)/V)/2 + a_1*(alpha_2 + w_g/V - (7*s*sigma/20 - s/2)*Derivative(q0(t), t)/V)*cos(alpha_2 - (7*s*sigma/20 - s/2)*Derivative(q0(t), t)/V))/20 + V**2*c*rho*s*sigma*(9*s*sigma/20 - s/2)*(C_Dmax*(1 - cos(2*alpha_2 + 2*w_g/V - 2*(9*s*sigma/20 - s/2)*Derivative(q0(t), t)/V))*sin(alpha_2 - (9*s*sigma/20 - s/2)*Derivative(q0(t), t)/V)/2 + a_1*(alpha_2 + w_g/V - (9*s*sigma/20 - s/2)*Derivative(q0(t), t)/V)*cos(alpha_2 - (9*s*sigma/20 - s/2)*Derivative(q0(t), t)/V))/20 - V*a_0*c*rho*s**3*(1 - sigma)**3*Derivative(q0(t), t)/24]])
	return e
