from sympy import *
def get_M():
	m_f = Symbol('m_f')
	s = Symbol('s')
	sigma = Symbol('sigma')
	t = Symbol('t')
	m_l = Symbol('m_l')
	m_w = Symbol('m_w')
	q1 = Function('q1')
	q2 = Function('q2')
	e = MutableDenseMatrix([[-m_f*s**2*sigma**2*cos(q1(t))/4 - m_f*s**2*sigma**2*cos(q2(t))/4 + 2*m_f*s**2*sigma**2/3 + m_f*s**2*sigma*cos(q1(t))/4 + m_f*s**2*sigma*cos(q2(t))/4 - m_f*s**2*sigma + m_f*s**2/2 + m_l*s**2*sigma**2/2 - m_l*s**2*sigma + m_l*s**2/2 + m_w*s**2*sigma**2/12 - m_w*s**2*sigma/6 + m_w*s**2/12, m_f*s**2*sigma**2*cos(q1(t))/8 - m_f*s**2*sigma**2/12 - m_f*s**2*sigma*cos(q1(t))/8, -m_f*s**2*sigma**2*cos(q2(t))/8 + m_f*s**2*sigma**2/12 + m_f*s**2*sigma*cos(q2(t))/8], [m_f*s**2*sigma**2*cos(q1(t))/8 - m_f*s**2*sigma**2/12 - m_f*s**2*sigma*cos(q1(t))/8, m_f*s**2*sigma**2/12, 0], [-m_f*s**2*sigma**2*cos(q2(t))/8 + m_f*s**2*sigma**2/12 + m_f*s**2*sigma*cos(q2(t))/8, 0, m_f*s**2*sigma**2/12]])
	return e
def get_f():
	g = Symbol('g')
	m_f = Symbol('m_f')
	s = Symbol('s')
	sigma = Symbol('sigma')
	t = Symbol('t')
	q0 = Function('q0')
	q1 = Function('q1')
	q2 = Function('q2')
	e = ImmutableDenseMatrix([[-g*m_f*(-s*sigma*(sin(q0(t))*sin(q1(t))/4 + cos(q0(t))*cos(q1(t))/4) - s*(Rational(1, 2) - sigma/2)*cos(q0(t))) - g*m_f*(-s*sigma*(sin(q0(t))*sin(q2(t))/4 - cos(q0(t))*cos(q2(t))/4) - s*(sigma/2 + Rational(-1, 2))*cos(q0(t))) + m_f*s**2*sigma**2*sin(q1(t))*Derivative(q0(t), t)*Derivative(q1(t), t)/4 - m_f*s**2*sigma**2*sin(q1(t))*Derivative(q1(t), t)**2/8 + m_f*s**2*sigma**2*sin(q2(t))*Derivative(q0(t), t)*Derivative(q2(t), t)/4 + m_f*s**2*sigma**2*sin(q2(t))*Derivative(q2(t), t)**2/8 - m_f*s**2*sigma*sin(q1(t))*Derivative(q0(t), t)*Derivative(q1(t), t)/4 + m_f*s**2*sigma*sin(q1(t))*Derivative(q1(t), t)**2/8 - m_f*s**2*sigma*sin(q2(t))*Derivative(q0(t), t)*Derivative(q2(t), t)/4 - m_f*s**2*sigma*sin(q2(t))*Derivative(q2(t), t)**2/8], [g*m_f*s*sigma*(-sin(q0(t))*sin(q1(t))/4 - cos(q0(t))*cos(q1(t))/4) - m_f*s**2*sigma**2*sin(q1(t))*Derivative(q0(t), t)*Derivative(q1(t), t)/8 + m_f*s**2*sigma*sin(q1(t))*Derivative(q0(t), t)*Derivative(q1(t), t)/8 - m_f*s**2*(3*sigma**2*sin(q1(t))*Derivative(q0(t), t)**2 - 3*sigma**2*sin(q1(t))*Derivative(q0(t), t)*Derivative(q1(t), t) - 3*sigma*sin(q1(t))*Derivative(q0(t), t)**2 + 3*sigma*sin(q1(t))*Derivative(q0(t), t)*Derivative(q1(t), t))/24], [g*m_f*s*sigma*(sin(q0(t))*sin(q2(t))/4 - cos(q0(t))*cos(q2(t))/4) + m_f*s**2*sigma**2*sin(q2(t))*Derivative(q0(t), t)*Derivative(q2(t), t)/8 - m_f*s**2*sigma*sin(q2(t))*Derivative(q0(t), t)*Derivative(q2(t), t)/8 - m_f*s**2*(3*sigma**2*sin(q2(t))*Derivative(q0(t), t)**2 + 3*sigma**2*sin(q2(t))*Derivative(q0(t), t)*Derivative(q2(t), t) - 3*sigma*sin(q2(t))*Derivative(q0(t), t)**2 - 3*sigma*sin(q2(t))*Derivative(q0(t), t)*Derivative(q2(t), t))/24]])
	return e
def get_T():
	m_f = Symbol('m_f')
	s = Symbol('s')
	sigma = Symbol('sigma')
	t = Symbol('t')
	m_l = Symbol('m_l')
	m_w = Symbol('m_w')
	q1 = Function('q1')
	q0 = Function('q0')
	q2 = Function('q2')
	e = m_f*s**2*(-3*sigma**2*cos(q1(t))*Derivative(q0(t), t)**2 + 3*sigma**2*cos(q1(t))*Derivative(q0(t), t)*Derivative(q1(t), t) + 4*sigma**2*Derivative(q0(t), t)**2 - 2*sigma**2*Derivative(q0(t), t)*Derivative(q1(t), t) + sigma**2*Derivative(q1(t), t)**2 + 3*sigma*cos(q1(t))*Derivative(q0(t), t)**2 - 3*sigma*cos(q1(t))*Derivative(q0(t), t)*Derivative(q1(t), t) - 6*sigma*Derivative(q0(t), t)**2 + 3*Derivative(q0(t), t)**2)/24 + m_f*s**2*(-3*sigma**2*cos(q2(t))*Derivative(q0(t), t)**2 - 3*sigma**2*cos(q2(t))*Derivative(q0(t), t)*Derivative(q2(t), t) + 4*sigma**2*Derivative(q0(t), t)**2 + 2*sigma**2*Derivative(q0(t), t)*Derivative(q2(t), t) + sigma**2*Derivative(q2(t), t)**2 + 3*sigma*cos(q2(t))*Derivative(q0(t), t)**2 + 3*sigma*cos(q2(t))*Derivative(q0(t), t)*Derivative(q2(t), t) - 6*sigma*Derivative(q0(t), t)**2 + 3*Derivative(q0(t), t)**2)/24 + m_l*s**2*(sigma - 1)**2*Derivative(q0(t), t)**2/4 + m_w*s**2*sigma**2*Derivative(q0(t), t)**2/24 - m_w*s**2*sigma*Derivative(q0(t), t)**2/12 + m_w*s**2*Derivative(q0(t), t)**2/24
	return e
def get_U():
	g = Symbol('g')
	m_f = Symbol('m_f')
	s = Symbol('s')
	sigma = Symbol('sigma')
	t = Symbol('t')
	q0 = Function('q0')
	q1 = Function('q1')
	q2 = Function('q2')
	e = g*m_f*(s*sigma*(sin(q0(t))*cos(q1(t)) - sin(q1(t))*cos(q0(t)))/4 + s*(1 - sigma)*sin(q0(t))/2) + g*m_f*(-s*sigma*(sin(q0(t))*cos(q2(t)) + sin(q2(t))*cos(q0(t)))/4 - s*(1 - sigma)*sin(q0(t))/2)
	return e
def get_Q():
	V = Symbol('V')
	c = Symbol('c')
	rho = Symbol('rho')
	s = Symbol('s')
	sigma = Symbol('sigma')
	eta_1 = Symbol('eta_1')
	C_Dmax = Symbol('C_Dmax')
	alpha_1 = Symbol('alpha_1')
	alpha_c = Symbol('alpha_c')
	eta_0 = Symbol('eta_0')
	w_g = Symbol('w_g')
	t = Symbol('t')
	a_0 = Symbol('a_0')
	alpha_2 = Symbol('alpha_2')
	q1 = Function('q1')
	q0 = Function('q0')
	q2 = Function('q2')
	e = ImmutableDenseMatrix([[V**2*c*rho*s*sigma*(eta_1/10 + Rational(9, 10))*(C_Dmax*(1 - cos(2*alpha_1 + 2*alpha_c + eta_0/5 + 2*w_g/V - 2*(-s*sigma*Derivative(q1(t), t)/20 + (-s*sigma*cos(q1(t))/2 + s*sigma/20 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_1 - (-s*sigma*Derivative(q1(t), t)/20 + (-s*sigma*cos(q1(t))/2 + s*sigma/20 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)/2 + a_0*(alpha_1 + alpha_c + eta_0/10 + w_g/V - (-s*sigma*Derivative(q1(t), t)/20 + (-s*sigma*cos(q1(t))/2 + s*sigma/20 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_1 - (-s*sigma*Derivative(q1(t), t)/20 + (-s*sigma*cos(q1(t))/2 + s*sigma/20 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*(-s*sigma*cos(q1(t))/2 + s*sigma/20 + s*cos(q1(t))/2)/20 + V**2*c*rho*s*sigma*(eta_1/10 + Rational(9, 10))*(C_Dmax*(1 - cos(2*alpha_2 + 2*alpha_c + eta_0/5 + 2*w_g/V - 2*(-s*sigma*Derivative(q2(t), t)/20 + (s*sigma*cos(q2(t))/2 - s*sigma/20 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_2 - (-s*sigma*Derivative(q2(t), t)/20 + (s*sigma*cos(q2(t))/2 - s*sigma/20 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)/2 + a_0*(alpha_2 + alpha_c + eta_0/10 + w_g/V - (-s*sigma*Derivative(q2(t), t)/20 + (s*sigma*cos(q2(t))/2 - s*sigma/20 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_2 - (-s*sigma*Derivative(q2(t), t)/20 + (s*sigma*cos(q2(t))/2 - s*sigma/20 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*(s*sigma*cos(q2(t))/2 - s*sigma/20 - s*cos(q2(t))/2)/20 + V**2*c*rho*s*sigma*(3*eta_1/10 + Rational(7, 10))*(C_Dmax*(1 - cos(2*alpha_1 + 2*alpha_c + 3*eta_0/5 + 2*w_g/V - 2*(-3*s*sigma*Derivative(q1(t), t)/20 + (-s*sigma*cos(q1(t))/2 + 3*s*sigma/20 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_1 - (-3*s*sigma*Derivative(q1(t), t)/20 + (-s*sigma*cos(q1(t))/2 + 3*s*sigma/20 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)/2 + a_0*(alpha_1 + alpha_c + 3*eta_0/10 + w_g/V - (-3*s*sigma*Derivative(q1(t), t)/20 + (-s*sigma*cos(q1(t))/2 + 3*s*sigma/20 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_1 - (-3*s*sigma*Derivative(q1(t), t)/20 + (-s*sigma*cos(q1(t))/2 + 3*s*sigma/20 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*(-s*sigma*cos(q1(t))/2 + 3*s*sigma/20 + s*cos(q1(t))/2)/20 + V**2*c*rho*s*sigma*(3*eta_1/10 + Rational(7, 10))*(C_Dmax*(1 - cos(2*alpha_2 + 2*alpha_c + 3*eta_0/5 + 2*w_g/V - 2*(-3*s*sigma*Derivative(q2(t), t)/20 + (s*sigma*cos(q2(t))/2 - 3*s*sigma/20 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_2 - (-3*s*sigma*Derivative(q2(t), t)/20 + (s*sigma*cos(q2(t))/2 - 3*s*sigma/20 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)/2 + a_0*(alpha_2 + alpha_c + 3*eta_0/10 + w_g/V - (-3*s*sigma*Derivative(q2(t), t)/20 + (s*sigma*cos(q2(t))/2 - 3*s*sigma/20 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_2 - (-3*s*sigma*Derivative(q2(t), t)/20 + (s*sigma*cos(q2(t))/2 - 3*s*sigma/20 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*(s*sigma*cos(q2(t))/2 - 3*s*sigma/20 - s*cos(q2(t))/2)/20 + V**2*c*rho*s*sigma*(eta_1/2 + Rational(1, 2))*(C_Dmax*(1 - cos(2*alpha_1 + 2*alpha_c + eta_0 + 2*w_g/V - 2*(-s*sigma*Derivative(q1(t), t)/4 + (-s*sigma*cos(q1(t))/2 + s*sigma/4 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_1 - (-s*sigma*Derivative(q1(t), t)/4 + (-s*sigma*cos(q1(t))/2 + s*sigma/4 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)/2 + a_0*(alpha_1 + alpha_c + eta_0/2 + w_g/V - (-s*sigma*Derivative(q1(t), t)/4 + (-s*sigma*cos(q1(t))/2 + s*sigma/4 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_1 - (-s*sigma*Derivative(q1(t), t)/4 + (-s*sigma*cos(q1(t))/2 + s*sigma/4 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*(-s*sigma*cos(q1(t))/2 + s*sigma/4 + s*cos(q1(t))/2)/20 + V**2*c*rho*s*sigma*(eta_1/2 + Rational(1, 2))*(C_Dmax*(1 - cos(2*alpha_2 + 2*alpha_c + eta_0 + 2*w_g/V - 2*(-s*sigma*Derivative(q2(t), t)/4 + (s*sigma*cos(q2(t))/2 - s*sigma/4 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_2 - (-s*sigma*Derivative(q2(t), t)/4 + (s*sigma*cos(q2(t))/2 - s*sigma/4 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)/2 + a_0*(alpha_2 + alpha_c + eta_0/2 + w_g/V - (-s*sigma*Derivative(q2(t), t)/4 + (s*sigma*cos(q2(t))/2 - s*sigma/4 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_2 - (-s*sigma*Derivative(q2(t), t)/4 + (s*sigma*cos(q2(t))/2 - s*sigma/4 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*(s*sigma*cos(q2(t))/2 - s*sigma/4 - s*cos(q2(t))/2)/20 + V**2*c*rho*s*sigma*(7*eta_1/10 + Rational(3, 10))*(C_Dmax*(1 - cos(2*alpha_1 + 2*alpha_c + 7*eta_0/5 + 2*w_g/V - 2*(-7*s*sigma*Derivative(q1(t), t)/20 + (-s*sigma*cos(q1(t))/2 + 7*s*sigma/20 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_1 - (-7*s*sigma*Derivative(q1(t), t)/20 + (-s*sigma*cos(q1(t))/2 + 7*s*sigma/20 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)/2 + a_0*(alpha_1 + alpha_c + 7*eta_0/10 + w_g/V - (-7*s*sigma*Derivative(q1(t), t)/20 + (-s*sigma*cos(q1(t))/2 + 7*s*sigma/20 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_1 - (-7*s*sigma*Derivative(q1(t), t)/20 + (-s*sigma*cos(q1(t))/2 + 7*s*sigma/20 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*(-s*sigma*cos(q1(t))/2 + 7*s*sigma/20 + s*cos(q1(t))/2)/20 + V**2*c*rho*s*sigma*(7*eta_1/10 + Rational(3, 10))*(C_Dmax*(1 - cos(2*alpha_2 + 2*alpha_c + 7*eta_0/5 + 2*w_g/V - 2*(-7*s*sigma*Derivative(q2(t), t)/20 + (s*sigma*cos(q2(t))/2 - 7*s*sigma/20 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_2 - (-7*s*sigma*Derivative(q2(t), t)/20 + (s*sigma*cos(q2(t))/2 - 7*s*sigma/20 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)/2 + a_0*(alpha_2 + alpha_c + 7*eta_0/10 + w_g/V - (-7*s*sigma*Derivative(q2(t), t)/20 + (s*sigma*cos(q2(t))/2 - 7*s*sigma/20 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_2 - (-7*s*sigma*Derivative(q2(t), t)/20 + (s*sigma*cos(q2(t))/2 - 7*s*sigma/20 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*(s*sigma*cos(q2(t))/2 - 7*s*sigma/20 - s*cos(q2(t))/2)/20 + V**2*c*rho*s*sigma*(9*eta_1/10 + Rational(1, 10))*(C_Dmax*(1 - cos(2*alpha_1 + 2*alpha_c + 9*eta_0/5 + 2*w_g/V - 2*(-9*s*sigma*Derivative(q1(t), t)/20 + (-s*sigma*cos(q1(t))/2 + 9*s*sigma/20 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_1 - (-9*s*sigma*Derivative(q1(t), t)/20 + (-s*sigma*cos(q1(t))/2 + 9*s*sigma/20 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)/2 + a_0*(alpha_1 + alpha_c + 9*eta_0/10 + w_g/V - (-9*s*sigma*Derivative(q1(t), t)/20 + (-s*sigma*cos(q1(t))/2 + 9*s*sigma/20 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_1 - (-9*s*sigma*Derivative(q1(t), t)/20 + (-s*sigma*cos(q1(t))/2 + 9*s*sigma/20 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*(-s*sigma*cos(q1(t))/2 + 9*s*sigma/20 + s*cos(q1(t))/2)/20 + V**2*c*rho*s*sigma*(9*eta_1/10 + Rational(1, 10))*(C_Dmax*(1 - cos(2*alpha_2 + 2*alpha_c + 9*eta_0/5 + 2*w_g/V - 2*(-9*s*sigma*Derivative(q2(t), t)/20 + (s*sigma*cos(q2(t))/2 - 9*s*sigma/20 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_2 - (-9*s*sigma*Derivative(q2(t), t)/20 + (s*sigma*cos(q2(t))/2 - 9*s*sigma/20 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)/2 + a_0*(alpha_2 + alpha_c + 9*eta_0/10 + w_g/V - (-9*s*sigma*Derivative(q2(t), t)/20 + (s*sigma*cos(q2(t))/2 - 9*s*sigma/20 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_2 - (-9*s*sigma*Derivative(q2(t), t)/20 + (s*sigma*cos(q2(t))/2 - 9*s*sigma/20 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*(s*sigma*cos(q2(t))/2 - 9*s*sigma/20 - s*cos(q2(t))/2)/20 - V*a_0*c*rho*s**3*(1 - sigma)**3*Derivative(q0(t), t)/24], [-V**2*c*rho*s**2*sigma**2*(eta_1/10 + Rational(9, 10))*(C_Dmax*(1 - cos(2*alpha_1 + 2*alpha_c + eta_0/5 + 2*w_g/V - 2*(-s*sigma*Derivative(q1(t), t)/20 + (-s*sigma*cos(q1(t))/2 + s*sigma/20 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_1 - (-s*sigma*Derivative(q1(t), t)/20 + (-s*sigma*cos(q1(t))/2 + s*sigma/20 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)/2 + a_0*(alpha_1 + alpha_c + eta_0/10 + w_g/V - (-s*sigma*Derivative(q1(t), t)/20 + (-s*sigma*cos(q1(t))/2 + s*sigma/20 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_1 - (-s*sigma*Derivative(q1(t), t)/20 + (-s*sigma*cos(q1(t))/2 + s*sigma/20 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))/400 - 3*V**2*c*rho*s**2*sigma**2*(3*eta_1/10 + Rational(7, 10))*(C_Dmax*(1 - cos(2*alpha_1 + 2*alpha_c + 3*eta_0/5 + 2*w_g/V - 2*(-3*s*sigma*Derivative(q1(t), t)/20 + (-s*sigma*cos(q1(t))/2 + 3*s*sigma/20 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_1 - (-3*s*sigma*Derivative(q1(t), t)/20 + (-s*sigma*cos(q1(t))/2 + 3*s*sigma/20 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)/2 + a_0*(alpha_1 + alpha_c + 3*eta_0/10 + w_g/V - (-3*s*sigma*Derivative(q1(t), t)/20 + (-s*sigma*cos(q1(t))/2 + 3*s*sigma/20 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_1 - (-3*s*sigma*Derivative(q1(t), t)/20 + (-s*sigma*cos(q1(t))/2 + 3*s*sigma/20 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))/400 - V**2*c*rho*s**2*sigma**2*(eta_1/2 + Rational(1, 2))*(C_Dmax*(1 - cos(2*alpha_1 + 2*alpha_c + eta_0 + 2*w_g/V - 2*(-s*sigma*Derivative(q1(t), t)/4 + (-s*sigma*cos(q1(t))/2 + s*sigma/4 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_1 - (-s*sigma*Derivative(q1(t), t)/4 + (-s*sigma*cos(q1(t))/2 + s*sigma/4 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)/2 + a_0*(alpha_1 + alpha_c + eta_0/2 + w_g/V - (-s*sigma*Derivative(q1(t), t)/4 + (-s*sigma*cos(q1(t))/2 + s*sigma/4 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_1 - (-s*sigma*Derivative(q1(t), t)/4 + (-s*sigma*cos(q1(t))/2 + s*sigma/4 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))/80 - 7*V**2*c*rho*s**2*sigma**2*(7*eta_1/10 + Rational(3, 10))*(C_Dmax*(1 - cos(2*alpha_1 + 2*alpha_c + 7*eta_0/5 + 2*w_g/V - 2*(-7*s*sigma*Derivative(q1(t), t)/20 + (-s*sigma*cos(q1(t))/2 + 7*s*sigma/20 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_1 - (-7*s*sigma*Derivative(q1(t), t)/20 + (-s*sigma*cos(q1(t))/2 + 7*s*sigma/20 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)/2 + a_0*(alpha_1 + alpha_c + 7*eta_0/10 + w_g/V - (-7*s*sigma*Derivative(q1(t), t)/20 + (-s*sigma*cos(q1(t))/2 + 7*s*sigma/20 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_1 - (-7*s*sigma*Derivative(q1(t), t)/20 + (-s*sigma*cos(q1(t))/2 + 7*s*sigma/20 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))/400 - 9*V**2*c*rho*s**2*sigma**2*(9*eta_1/10 + Rational(1, 10))*(C_Dmax*(1 - cos(2*alpha_1 + 2*alpha_c + 9*eta_0/5 + 2*w_g/V - 2*(-9*s*sigma*Derivative(q1(t), t)/20 + (-s*sigma*cos(q1(t))/2 + 9*s*sigma/20 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_1 - (-9*s*sigma*Derivative(q1(t), t)/20 + (-s*sigma*cos(q1(t))/2 + 9*s*sigma/20 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)/2 + a_0*(alpha_1 + alpha_c + 9*eta_0/10 + w_g/V - (-9*s*sigma*Derivative(q1(t), t)/20 + (-s*sigma*cos(q1(t))/2 + 9*s*sigma/20 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_1 - (-9*s*sigma*Derivative(q1(t), t)/20 + (-s*sigma*cos(q1(t))/2 + 9*s*sigma/20 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))/400], [-V**2*c*rho*s**2*sigma**2*(eta_1/10 + Rational(9, 10))*(C_Dmax*(1 - cos(2*alpha_2 + 2*alpha_c + eta_0/5 + 2*w_g/V - 2*(-s*sigma*Derivative(q2(t), t)/20 + (s*sigma*cos(q2(t))/2 - s*sigma/20 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_2 - (-s*sigma*Derivative(q2(t), t)/20 + (s*sigma*cos(q2(t))/2 - s*sigma/20 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)/2 + a_0*(alpha_2 + alpha_c + eta_0/10 + w_g/V - (-s*sigma*Derivative(q2(t), t)/20 + (s*sigma*cos(q2(t))/2 - s*sigma/20 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_2 - (-s*sigma*Derivative(q2(t), t)/20 + (s*sigma*cos(q2(t))/2 - s*sigma/20 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))/400 - 3*V**2*c*rho*s**2*sigma**2*(3*eta_1/10 + Rational(7, 10))*(C_Dmax*(1 - cos(2*alpha_2 + 2*alpha_c + 3*eta_0/5 + 2*w_g/V - 2*(-3*s*sigma*Derivative(q2(t), t)/20 + (s*sigma*cos(q2(t))/2 - 3*s*sigma/20 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_2 - (-3*s*sigma*Derivative(q2(t), t)/20 + (s*sigma*cos(q2(t))/2 - 3*s*sigma/20 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)/2 + a_0*(alpha_2 + alpha_c + 3*eta_0/10 + w_g/V - (-3*s*sigma*Derivative(q2(t), t)/20 + (s*sigma*cos(q2(t))/2 - 3*s*sigma/20 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_2 - (-3*s*sigma*Derivative(q2(t), t)/20 + (s*sigma*cos(q2(t))/2 - 3*s*sigma/20 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))/400 - V**2*c*rho*s**2*sigma**2*(eta_1/2 + Rational(1, 2))*(C_Dmax*(1 - cos(2*alpha_2 + 2*alpha_c + eta_0 + 2*w_g/V - 2*(-s*sigma*Derivative(q2(t), t)/4 + (s*sigma*cos(q2(t))/2 - s*sigma/4 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_2 - (-s*sigma*Derivative(q2(t), t)/4 + (s*sigma*cos(q2(t))/2 - s*sigma/4 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)/2 + a_0*(alpha_2 + alpha_c + eta_0/2 + w_g/V - (-s*sigma*Derivative(q2(t), t)/4 + (s*sigma*cos(q2(t))/2 - s*sigma/4 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_2 - (-s*sigma*Derivative(q2(t), t)/4 + (s*sigma*cos(q2(t))/2 - s*sigma/4 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))/80 - 7*V**2*c*rho*s**2*sigma**2*(7*eta_1/10 + Rational(3, 10))*(C_Dmax*(1 - cos(2*alpha_2 + 2*alpha_c + 7*eta_0/5 + 2*w_g/V - 2*(-7*s*sigma*Derivative(q2(t), t)/20 + (s*sigma*cos(q2(t))/2 - 7*s*sigma/20 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_2 - (-7*s*sigma*Derivative(q2(t), t)/20 + (s*sigma*cos(q2(t))/2 - 7*s*sigma/20 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)/2 + a_0*(alpha_2 + alpha_c + 7*eta_0/10 + w_g/V - (-7*s*sigma*Derivative(q2(t), t)/20 + (s*sigma*cos(q2(t))/2 - 7*s*sigma/20 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_2 - (-7*s*sigma*Derivative(q2(t), t)/20 + (s*sigma*cos(q2(t))/2 - 7*s*sigma/20 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))/400 - 9*V**2*c*rho*s**2*sigma**2*(9*eta_1/10 + Rational(1, 10))*(C_Dmax*(1 - cos(2*alpha_2 + 2*alpha_c + 9*eta_0/5 + 2*w_g/V - 2*(-9*s*sigma*Derivative(q2(t), t)/20 + (s*sigma*cos(q2(t))/2 - 9*s*sigma/20 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_2 - (-9*s*sigma*Derivative(q2(t), t)/20 + (s*sigma*cos(q2(t))/2 - 9*s*sigma/20 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)/2 + a_0*(alpha_2 + alpha_c + 9*eta_0/10 + w_g/V - (-9*s*sigma*Derivative(q2(t), t)/20 + (s*sigma*cos(q2(t))/2 - 9*s*sigma/20 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_2 - (-9*s*sigma*Derivative(q2(t), t)/20 + (s*sigma*cos(q2(t))/2 - 9*s*sigma/20 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))/400]])
	return e
