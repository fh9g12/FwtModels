from sympy import *
def get_M():
	I_xxf = Symbol('I_xxf')
	I_xxw = Symbol('I_xxw')
	eta_0 = Symbol('eta_0')
	m_w = Symbol('m_w')
	eta_1 = Symbol('eta_1')
	eta_2 = Symbol('eta_2')
	m_f = Symbol('m_f')
	s = Symbol('s')
	sigma = Symbol('sigma')
	t = Symbol('t')
	q1 = Function('q1')
	q2 = Function('q2')
	e = MutableDenseMatrix([[2*I_xxf + I_xxw + eta_0**2*m_w + eta_1**2*m_w + 2*eta_2**2*m_f - eta_2*m_f*s*sigma*cos(q1(t)) - eta_2*m_f*s*sigma*cos(q2(t)) + eta_2*m_f*s*cos(q1(t)) + eta_2*m_f*s*cos(q2(t)) + m_f*s**2*sigma**2/2 - m_f*s**2*sigma + m_f*s**2/2, -I_xxf - eta_2**2*m_f + eta_2*m_f*s*sigma*cos(q1(t))/2 - eta_2*m_f*s*cos(q1(t))/2, I_xxf + eta_2**2*m_f - eta_2*m_f*s*sigma*cos(q2(t))/2 + eta_2*m_f*s*cos(q2(t))/2], [-I_xxf - eta_2**2*m_f + eta_2*m_f*s*sigma*cos(q1(t))/2 - eta_2*m_f*s*cos(q1(t))/2, I_xxf + eta_2**2*m_f, 0], [I_xxf + eta_2**2*m_f - eta_2*m_f*s*sigma*cos(q2(t))/2 + eta_2*m_f*s*cos(q2(t))/2, 0, I_xxf + eta_2**2*m_f]])
	return e
def get_f():
	eta_2 = Symbol('eta_2')
	m_f = Symbol('m_f')
	s = Symbol('s')
	sigma = Symbol('sigma')
	t = Symbol('t')
	g = Symbol('g')
	m_w = Symbol('m_w')
	eta_0 = Symbol('eta_0')
	eta_1 = Symbol('eta_1')
	q1 = Function('q1')
	q0 = Function('q0')
	q2 = Function('q2')
	e = ImmutableDenseMatrix([[eta_2*m_f*s*sigma*sin(q1(t))*Derivative(q0(t), t)*Derivative(q1(t), t) - eta_2*m_f*s*sigma*sin(q1(t))*Derivative(q1(t), t)**2/2 + eta_2*m_f*s*sigma*sin(q2(t))*Derivative(q0(t), t)*Derivative(q2(t), t) + eta_2*m_f*s*sigma*sin(q2(t))*Derivative(q2(t), t)**2/2 - eta_2*m_f*s*sin(q1(t))*Derivative(q0(t), t)*Derivative(q1(t), t) + eta_2*m_f*s*sin(q1(t))*Derivative(q1(t), t)**2/2 - eta_2*m_f*s*sin(q2(t))*Derivative(q0(t), t)*Derivative(q2(t), t) - eta_2*m_f*s*sin(q2(t))*Derivative(q2(t), t)**2/2 - g*m_f*(-eta_2*(sin(q0(t))*sin(q1(t)) + cos(q0(t))*cos(q1(t))) - s*(Rational(1, 2) - sigma/2)*cos(q0(t))) - g*m_f*(-eta_2*(sin(q0(t))*sin(q2(t)) - cos(q0(t))*cos(q2(t))) - s*(sigma/2 + Rational(-1, 2))*cos(q0(t))) - g*m_w*(-eta_0*cos(q0(t)) + eta_1*sin(q0(t)))], [eta_2*g*m_f*(-sin(q0(t))*sin(q1(t)) - cos(q0(t))*cos(q1(t))) - eta_2*m_f*s*sigma*sin(q1(t))*Derivative(q0(t), t)**2/2 + eta_2*m_f*s*sin(q1(t))*Derivative(q0(t), t)**2/2], [eta_2*g*m_f*(sin(q0(t))*sin(q2(t)) - cos(q0(t))*cos(q2(t))) - eta_2*m_f*s*sigma*sin(q2(t))*Derivative(q0(t), t)**2/2 + eta_2*m_f*s*sin(q2(t))*Derivative(q0(t), t)**2/2]])
	return e
def get_T():
	I_xxf = Symbol('I_xxf')
	t = Symbol('t')
	I_xxw = Symbol('I_xxw')
	eta_0 = Symbol('eta_0')
	m_w = Symbol('m_w')
	eta_1 = Symbol('eta_1')
	eta_2 = Symbol('eta_2')
	m_f = Symbol('m_f')
	s = Symbol('s')
	sigma = Symbol('sigma')
	q0 = Function('q0')
	q1 = Function('q1')
	q2 = Function('q2')
	e = I_xxf*Derivative(q0(t), t)**2 - I_xxf*Derivative(q0(t), t)*Derivative(q1(t), t) + I_xxf*Derivative(q0(t), t)*Derivative(q2(t), t) + I_xxf*Derivative(q1(t), t)**2/2 + I_xxf*Derivative(q2(t), t)**2/2 + I_xxw*Derivative(q0(t), t)**2/2 + eta_0**2*m_w*Derivative(q0(t), t)**2/2 + eta_1**2*m_w*Derivative(q0(t), t)**2/2 + eta_2**2*m_f*Derivative(q0(t), t)**2 - eta_2**2*m_f*Derivative(q0(t), t)*Derivative(q1(t), t) + eta_2**2*m_f*Derivative(q0(t), t)*Derivative(q2(t), t) + eta_2**2*m_f*Derivative(q1(t), t)**2/2 + eta_2**2*m_f*Derivative(q2(t), t)**2/2 - eta_2*m_f*s*sigma*cos(q1(t))*Derivative(q0(t), t)**2/2 + eta_2*m_f*s*sigma*cos(q1(t))*Derivative(q0(t), t)*Derivative(q1(t), t)/2 - eta_2*m_f*s*sigma*cos(q2(t))*Derivative(q0(t), t)**2/2 - eta_2*m_f*s*sigma*cos(q2(t))*Derivative(q0(t), t)*Derivative(q2(t), t)/2 + eta_2*m_f*s*cos(q1(t))*Derivative(q0(t), t)**2/2 - eta_2*m_f*s*cos(q1(t))*Derivative(q0(t), t)*Derivative(q1(t), t)/2 + eta_2*m_f*s*cos(q2(t))*Derivative(q0(t), t)**2/2 + eta_2*m_f*s*cos(q2(t))*Derivative(q0(t), t)*Derivative(q2(t), t)/2 + m_f*s**2*sigma**2*Derivative(q0(t), t)**2/4 - m_f*s**2*sigma*Derivative(q0(t), t)**2/2 + m_f*s**2*Derivative(q0(t), t)**2/4
	return e
def get_U():
	g = Symbol('g')
	m_f = Symbol('m_f')
	eta_2 = Symbol('eta_2')
	t = Symbol('t')
	s = Symbol('s')
	sigma = Symbol('sigma')
	m_w = Symbol('m_w')
	eta_0 = Symbol('eta_0')
	eta_1 = Symbol('eta_1')
	q0 = Function('q0')
	q1 = Function('q1')
	q2 = Function('q2')
	e = g*m_f*(eta_2*(sin(q0(t))*cos(q1(t)) - sin(q1(t))*cos(q0(t))) + s*(1 - sigma)*sin(q0(t))/2) + g*m_f*(-eta_2*(sin(q0(t))*cos(q2(t)) + sin(q2(t))*cos(q0(t))) - s*(1 - sigma)*sin(q0(t))/2) + g*m_w*(eta_0*sin(q0(t)) + eta_1*cos(q0(t)))
	return e
def get_Q():
	V = Symbol('V')
	c = Symbol('c')
	rho = Symbol('rho')
	s = Symbol('s')
	sigma = Symbol('sigma')
	C_Dmax = Symbol('C_Dmax')
	alpha_1 = Symbol('alpha_1')
	alpha_c = Symbol('alpha_c')
	w_g = Symbol('w_g')
	Lambda = Symbol('Lambda')
	t = Symbol('t')
	a29 = Symbol('a29')
	a28 = Symbol('a28')
	a27 = Symbol('a27')
	a26 = Symbol('a26')
	a25 = Symbol('a25')
	a24 = Symbol('a24')
	a23 = Symbol('a23')
	a22 = Symbol('a22')
	a21 = Symbol('a21')
	a20 = Symbol('a20')
	alpha_2 = Symbol('alpha_2')
	a9 = Symbol('a9')
	a8 = Symbol('a8')
	a7 = Symbol('a7')
	a6 = Symbol('a6')
	a5 = Symbol('a5')
	a4 = Symbol('a4')
	a3 = Symbol('a3')
	a2 = Symbol('a2')
	a1 = Symbol('a1')
	a0 = Symbol('a0')
	a10 = Symbol('a10')
	a11 = Symbol('a11')
	a12 = Symbol('a12')
	a13 = Symbol('a13')
	a14 = Symbol('a14')
	a15 = Symbol('a15')
	a16 = Symbol('a16')
	a17 = Symbol('a17')
	a18 = Symbol('a18')
	a19 = Symbol('a19')
	q1 = Function('q1')
	q0 = Function('q0')
	q2 = Function('q2')
	e = ImmutableDenseMatrix([[V**2*c*rho*s*sigma*(C_Dmax*(1 - cos(2*alpha_1 + 2*alpha_c + 2*w_g/V - 2*(-19*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (19*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_1 - (-19*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (19*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)/2 + a29*(alpha_1 + alpha_c + w_g/V - (-19*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (19*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_1 - (-19*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (19*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*(19*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)/40 + V**2*c*rho*s*sigma*(C_Dmax*(1 - cos(2*alpha_1 + 2*alpha_c + 2*w_g/V - 2*(-17*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (17*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_1 - (-17*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (17*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)/2 + a28*(alpha_1 + alpha_c + w_g/V - (-17*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (17*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_1 - (-17*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (17*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*(17*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)/40 + V**2*c*rho*s*sigma*(C_Dmax*(1 - cos(2*alpha_1 + 2*alpha_c + 2*w_g/V - 2*(-3*s*sigma*cos(Lambda)*Derivative(q1(t), t)/8 + (3*s*sigma*cos(Lambda)/8 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_1 - (-3*s*sigma*cos(Lambda)*Derivative(q1(t), t)/8 + (3*s*sigma*cos(Lambda)/8 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)/2 + a27*(alpha_1 + alpha_c + w_g/V - (-3*s*sigma*cos(Lambda)*Derivative(q1(t), t)/8 + (3*s*sigma*cos(Lambda)/8 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_1 - (-3*s*sigma*cos(Lambda)*Derivative(q1(t), t)/8 + (3*s*sigma*cos(Lambda)/8 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*(3*s*sigma*cos(Lambda)/8 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)/40 + V**2*c*rho*s*sigma*(C_Dmax*(1 - cos(2*alpha_1 + 2*alpha_c + 2*w_g/V - 2*(-13*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (13*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_1 - (-13*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (13*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)/2 + a26*(alpha_1 + alpha_c + w_g/V - (-13*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (13*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_1 - (-13*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (13*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*(13*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)/40 + V**2*c*rho*s*sigma*(C_Dmax*(1 - cos(2*alpha_1 + 2*alpha_c + 2*w_g/V - 2*(-11*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (11*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_1 - (-11*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (11*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)/2 + a25*(alpha_1 + alpha_c + w_g/V - (-11*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (11*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_1 - (-11*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (11*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*(11*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)/40 + V**2*c*rho*s*sigma*(C_Dmax*(1 - cos(2*alpha_1 + 2*alpha_c + 2*w_g/V - 2*(-9*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (9*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_1 - (-9*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (9*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)/2 + a24*(alpha_1 + alpha_c + w_g/V - (-9*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (9*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_1 - (-9*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (9*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*(9*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)/40 + V**2*c*rho*s*sigma*(C_Dmax*(1 - cos(2*alpha_1 + 2*alpha_c + 2*w_g/V - 2*(-7*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (7*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_1 - (-7*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (7*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)/2 + a23*(alpha_1 + alpha_c + w_g/V - (-7*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (7*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_1 - (-7*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (7*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*(7*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)/40 + V**2*c*rho*s*sigma*(C_Dmax*(1 - cos(2*alpha_1 + 2*alpha_c + 2*w_g/V - 2*(-s*sigma*cos(Lambda)*Derivative(q1(t), t)/8 + (s*sigma*cos(Lambda)/8 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_1 - (-s*sigma*cos(Lambda)*Derivative(q1(t), t)/8 + (s*sigma*cos(Lambda)/8 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)/2 + a22*(alpha_1 + alpha_c + w_g/V - (-s*sigma*cos(Lambda)*Derivative(q1(t), t)/8 + (s*sigma*cos(Lambda)/8 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_1 - (-s*sigma*cos(Lambda)*Derivative(q1(t), t)/8 + (s*sigma*cos(Lambda)/8 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*(s*sigma*cos(Lambda)/8 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)/40 + V**2*c*rho*s*sigma*(C_Dmax*(1 - cos(2*alpha_1 + 2*alpha_c + 2*w_g/V - 2*(-3*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (3*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_1 - (-3*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (3*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)/2 + a21*(alpha_1 + alpha_c + w_g/V - (-3*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (3*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_1 - (-3*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (3*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*(3*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)/40 + V**2*c*rho*s*sigma*(C_Dmax*(1 - cos(2*alpha_1 + 2*alpha_c + 2*w_g/V - 2*(-s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_1 - (-s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)/2 + a20*(alpha_1 + alpha_c + w_g/V - (-s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_1 - (-s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*(s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)/40 + V**2*c*rho*s*sigma*(C_Dmax*(1 - cos(2*alpha_2 + 2*alpha_c + 2*w_g/V - 2*(-19*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-19*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_2 - (-19*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-19*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)/2 + a9*(alpha_2 + alpha_c + w_g/V - (-19*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-19*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_2 - (-19*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-19*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*(-19*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)/40 + V**2*c*rho*s*sigma*(C_Dmax*(1 - cos(2*alpha_2 + 2*alpha_c + 2*w_g/V - 2*(-17*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-17*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_2 - (-17*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-17*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)/2 + a8*(alpha_2 + alpha_c + w_g/V - (-17*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-17*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_2 - (-17*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-17*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*(-17*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)/40 + V**2*c*rho*s*sigma*(C_Dmax*(1 - cos(2*alpha_2 + 2*alpha_c + 2*w_g/V - 2*(-3*s*sigma*cos(Lambda)*Derivative(q2(t), t)/8 + (-3*s*sigma*cos(Lambda)/8 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_2 - (-3*s*sigma*cos(Lambda)*Derivative(q2(t), t)/8 + (-3*s*sigma*cos(Lambda)/8 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)/2 + a7*(alpha_2 + alpha_c + w_g/V - (-3*s*sigma*cos(Lambda)*Derivative(q2(t), t)/8 + (-3*s*sigma*cos(Lambda)/8 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_2 - (-3*s*sigma*cos(Lambda)*Derivative(q2(t), t)/8 + (-3*s*sigma*cos(Lambda)/8 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*(-3*s*sigma*cos(Lambda)/8 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)/40 + V**2*c*rho*s*sigma*(C_Dmax*(1 - cos(2*alpha_2 + 2*alpha_c + 2*w_g/V - 2*(-13*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-13*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_2 - (-13*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-13*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)/2 + a6*(alpha_2 + alpha_c + w_g/V - (-13*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-13*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_2 - (-13*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-13*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*(-13*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)/40 + V**2*c*rho*s*sigma*(C_Dmax*(1 - cos(2*alpha_2 + 2*alpha_c + 2*w_g/V - 2*(-11*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-11*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_2 - (-11*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-11*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)/2 + a5*(alpha_2 + alpha_c + w_g/V - (-11*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-11*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_2 - (-11*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-11*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*(-11*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)/40 + V**2*c*rho*s*sigma*(C_Dmax*(1 - cos(2*alpha_2 + 2*alpha_c + 2*w_g/V - 2*(-9*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-9*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_2 - (-9*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-9*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)/2 + a4*(alpha_2 + alpha_c + w_g/V - (-9*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-9*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_2 - (-9*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-9*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*(-9*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)/40 + V**2*c*rho*s*sigma*(C_Dmax*(1 - cos(2*alpha_2 + 2*alpha_c + 2*w_g/V - 2*(-7*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-7*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_2 - (-7*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-7*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)/2 + a3*(alpha_2 + alpha_c + w_g/V - (-7*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-7*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_2 - (-7*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-7*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*(-7*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)/40 + V**2*c*rho*s*sigma*(C_Dmax*(1 - cos(2*alpha_2 + 2*alpha_c + 2*w_g/V - 2*(-s*sigma*cos(Lambda)*Derivative(q2(t), t)/8 + (-s*sigma*cos(Lambda)/8 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_2 - (-s*sigma*cos(Lambda)*Derivative(q2(t), t)/8 + (-s*sigma*cos(Lambda)/8 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)/2 + a2*(alpha_2 + alpha_c + w_g/V - (-s*sigma*cos(Lambda)*Derivative(q2(t), t)/8 + (-s*sigma*cos(Lambda)/8 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_2 - (-s*sigma*cos(Lambda)*Derivative(q2(t), t)/8 + (-s*sigma*cos(Lambda)/8 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*(-s*sigma*cos(Lambda)/8 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)/40 + V**2*c*rho*s*sigma*(C_Dmax*(1 - cos(2*alpha_2 + 2*alpha_c + 2*w_g/V - 2*(-3*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-3*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_2 - (-3*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-3*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)/2 + a1*(alpha_2 + alpha_c + w_g/V - (-3*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-3*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_2 - (-3*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-3*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*(-3*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)/40 + V**2*c*rho*s*sigma*(C_Dmax*(1 - cos(2*alpha_2 + 2*alpha_c + 2*w_g/V - 2*(-s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_2 - (-s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)/2 + a0*(alpha_2 + alpha_c + w_g/V - (-s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_2 - (-s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*(-s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)/40 - 81*V*c*rho*s**3*a10*(1 - sigma)**3*Derivative(q0(t), t)/8000 - 49*V*c*rho*s**3*a11*(1 - sigma)**3*Derivative(q0(t), t)/8000 - V*c*rho*s**3*a12*(1 - sigma)**3*Derivative(q0(t), t)/320 - 9*V*c*rho*s**3*a13*(1 - sigma)**3*Derivative(q0(t), t)/8000 - V*c*rho*s**3*a14*(1 - sigma)**3*Derivative(q0(t), t)/8000 - V*c*rho*s**3*a15*(1 - sigma)**3*Derivative(q0(t), t)/8000 - 9*V*c*rho*s**3*a16*(1 - sigma)**3*Derivative(q0(t), t)/8000 - V*c*rho*s**3*a17*(1 - sigma)**3*Derivative(q0(t), t)/320 - 49*V*c*rho*s**3*a18*(1 - sigma)**3*Derivative(q0(t), t)/8000 - 81*V*c*rho*s**3*a19*(1 - sigma)**3*Derivative(q0(t), t)/8000], [-19*V**2*c*rho*s**2*sigma**2*(C_Dmax*(1 - cos(2*alpha_1 + 2*alpha_c + 2*w_g/V - 2*(-19*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (19*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_1 - (-19*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (19*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)/2 + a29*(alpha_1 + alpha_c + w_g/V - (-19*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (19*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_1 - (-19*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (19*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*cos(Lambda)/1600 - 17*V**2*c*rho*s**2*sigma**2*(C_Dmax*(1 - cos(2*alpha_1 + 2*alpha_c + 2*w_g/V - 2*(-17*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (17*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_1 - (-17*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (17*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)/2 + a28*(alpha_1 + alpha_c + w_g/V - (-17*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (17*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_1 - (-17*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (17*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*cos(Lambda)/1600 - 3*V**2*c*rho*s**2*sigma**2*(C_Dmax*(1 - cos(2*alpha_1 + 2*alpha_c + 2*w_g/V - 2*(-3*s*sigma*cos(Lambda)*Derivative(q1(t), t)/8 + (3*s*sigma*cos(Lambda)/8 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_1 - (-3*s*sigma*cos(Lambda)*Derivative(q1(t), t)/8 + (3*s*sigma*cos(Lambda)/8 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)/2 + a27*(alpha_1 + alpha_c + w_g/V - (-3*s*sigma*cos(Lambda)*Derivative(q1(t), t)/8 + (3*s*sigma*cos(Lambda)/8 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_1 - (-3*s*sigma*cos(Lambda)*Derivative(q1(t), t)/8 + (3*s*sigma*cos(Lambda)/8 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*cos(Lambda)/320 - 13*V**2*c*rho*s**2*sigma**2*(C_Dmax*(1 - cos(2*alpha_1 + 2*alpha_c + 2*w_g/V - 2*(-13*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (13*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_1 - (-13*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (13*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)/2 + a26*(alpha_1 + alpha_c + w_g/V - (-13*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (13*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_1 - (-13*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (13*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*cos(Lambda)/1600 - 11*V**2*c*rho*s**2*sigma**2*(C_Dmax*(1 - cos(2*alpha_1 + 2*alpha_c + 2*w_g/V - 2*(-11*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (11*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_1 - (-11*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (11*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)/2 + a25*(alpha_1 + alpha_c + w_g/V - (-11*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (11*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_1 - (-11*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (11*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*cos(Lambda)/1600 - 9*V**2*c*rho*s**2*sigma**2*(C_Dmax*(1 - cos(2*alpha_1 + 2*alpha_c + 2*w_g/V - 2*(-9*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (9*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_1 - (-9*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (9*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)/2 + a24*(alpha_1 + alpha_c + w_g/V - (-9*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (9*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_1 - (-9*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (9*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*cos(Lambda)/1600 - 7*V**2*c*rho*s**2*sigma**2*(C_Dmax*(1 - cos(2*alpha_1 + 2*alpha_c + 2*w_g/V - 2*(-7*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (7*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_1 - (-7*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (7*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)/2 + a23*(alpha_1 + alpha_c + w_g/V - (-7*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (7*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_1 - (-7*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (7*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*cos(Lambda)/1600 - V**2*c*rho*s**2*sigma**2*(C_Dmax*(1 - cos(2*alpha_1 + 2*alpha_c + 2*w_g/V - 2*(-s*sigma*cos(Lambda)*Derivative(q1(t), t)/8 + (s*sigma*cos(Lambda)/8 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_1 - (-s*sigma*cos(Lambda)*Derivative(q1(t), t)/8 + (s*sigma*cos(Lambda)/8 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)/2 + a22*(alpha_1 + alpha_c + w_g/V - (-s*sigma*cos(Lambda)*Derivative(q1(t), t)/8 + (s*sigma*cos(Lambda)/8 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_1 - (-s*sigma*cos(Lambda)*Derivative(q1(t), t)/8 + (s*sigma*cos(Lambda)/8 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*cos(Lambda)/320 - 3*V**2*c*rho*s**2*sigma**2*(C_Dmax*(1 - cos(2*alpha_1 + 2*alpha_c + 2*w_g/V - 2*(-3*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (3*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_1 - (-3*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (3*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)/2 + a21*(alpha_1 + alpha_c + w_g/V - (-3*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (3*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_1 - (-3*s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (3*s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*cos(Lambda)/1600 - V**2*c*rho*s**2*sigma**2*(C_Dmax*(1 - cos(2*alpha_1 + 2*alpha_c + 2*w_g/V - 2*(-s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_1 - (-s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)/2 + a20*(alpha_1 + alpha_c + w_g/V - (-s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_1 - (-s*sigma*cos(Lambda)*Derivative(q1(t), t)/40 + (s*sigma*cos(Lambda)/40 - s*sigma*cos(q1(t))/2 + s*cos(q1(t))/2)*Derivative(q0(t), t))/V))*cos(Lambda)/1600], [-19*V**2*c*rho*s**2*sigma**2*(C_Dmax*(1 - cos(2*alpha_2 + 2*alpha_c + 2*w_g/V - 2*(-19*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-19*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_2 - (-19*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-19*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)/2 + a9*(alpha_2 + alpha_c + w_g/V - (-19*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-19*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_2 - (-19*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-19*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*cos(Lambda)/1600 - 17*V**2*c*rho*s**2*sigma**2*(C_Dmax*(1 - cos(2*alpha_2 + 2*alpha_c + 2*w_g/V - 2*(-17*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-17*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_2 - (-17*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-17*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)/2 + a8*(alpha_2 + alpha_c + w_g/V - (-17*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-17*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_2 - (-17*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-17*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*cos(Lambda)/1600 - 3*V**2*c*rho*s**2*sigma**2*(C_Dmax*(1 - cos(2*alpha_2 + 2*alpha_c + 2*w_g/V - 2*(-3*s*sigma*cos(Lambda)*Derivative(q2(t), t)/8 + (-3*s*sigma*cos(Lambda)/8 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_2 - (-3*s*sigma*cos(Lambda)*Derivative(q2(t), t)/8 + (-3*s*sigma*cos(Lambda)/8 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)/2 + a7*(alpha_2 + alpha_c + w_g/V - (-3*s*sigma*cos(Lambda)*Derivative(q2(t), t)/8 + (-3*s*sigma*cos(Lambda)/8 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_2 - (-3*s*sigma*cos(Lambda)*Derivative(q2(t), t)/8 + (-3*s*sigma*cos(Lambda)/8 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*cos(Lambda)/320 - 13*V**2*c*rho*s**2*sigma**2*(C_Dmax*(1 - cos(2*alpha_2 + 2*alpha_c + 2*w_g/V - 2*(-13*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-13*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_2 - (-13*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-13*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)/2 + a6*(alpha_2 + alpha_c + w_g/V - (-13*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-13*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_2 - (-13*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-13*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*cos(Lambda)/1600 - 11*V**2*c*rho*s**2*sigma**2*(C_Dmax*(1 - cos(2*alpha_2 + 2*alpha_c + 2*w_g/V - 2*(-11*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-11*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_2 - (-11*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-11*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)/2 + a5*(alpha_2 + alpha_c + w_g/V - (-11*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-11*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_2 - (-11*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-11*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*cos(Lambda)/1600 - 9*V**2*c*rho*s**2*sigma**2*(C_Dmax*(1 - cos(2*alpha_2 + 2*alpha_c + 2*w_g/V - 2*(-9*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-9*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_2 - (-9*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-9*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)/2 + a4*(alpha_2 + alpha_c + w_g/V - (-9*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-9*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_2 - (-9*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-9*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*cos(Lambda)/1600 - 7*V**2*c*rho*s**2*sigma**2*(C_Dmax*(1 - cos(2*alpha_2 + 2*alpha_c + 2*w_g/V - 2*(-7*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-7*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_2 - (-7*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-7*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)/2 + a3*(alpha_2 + alpha_c + w_g/V - (-7*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-7*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_2 - (-7*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-7*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*cos(Lambda)/1600 - V**2*c*rho*s**2*sigma**2*(C_Dmax*(1 - cos(2*alpha_2 + 2*alpha_c + 2*w_g/V - 2*(-s*sigma*cos(Lambda)*Derivative(q2(t), t)/8 + (-s*sigma*cos(Lambda)/8 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_2 - (-s*sigma*cos(Lambda)*Derivative(q2(t), t)/8 + (-s*sigma*cos(Lambda)/8 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)/2 + a2*(alpha_2 + alpha_c + w_g/V - (-s*sigma*cos(Lambda)*Derivative(q2(t), t)/8 + (-s*sigma*cos(Lambda)/8 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_2 - (-s*sigma*cos(Lambda)*Derivative(q2(t), t)/8 + (-s*sigma*cos(Lambda)/8 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*cos(Lambda)/320 - 3*V**2*c*rho*s**2*sigma**2*(C_Dmax*(1 - cos(2*alpha_2 + 2*alpha_c + 2*w_g/V - 2*(-3*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-3*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_2 - (-3*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-3*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)/2 + a1*(alpha_2 + alpha_c + w_g/V - (-3*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-3*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_2 - (-3*s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-3*s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*cos(Lambda)/1600 - V**2*c*rho*s**2*sigma**2*(C_Dmax*(1 - cos(2*alpha_2 + 2*alpha_c + 2*w_g/V - 2*(-s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*sin(alpha_2 - (-s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)/2 + a0*(alpha_2 + alpha_c + w_g/V - (-s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V)*cos(alpha_2 - (-s*sigma*cos(Lambda)*Derivative(q2(t), t)/40 + (-s*sigma*cos(Lambda)/40 + s*sigma*cos(q2(t))/2 - s*cos(q2(t))/2)*Derivative(q0(t), t))/V))*cos(Lambda)/1600]])
	return e
