from sympy import *
def get_M():
	m_f = Symbol('m_f')
	s = Symbol('s')
	sigma = Symbol('sigma')
	t = Symbol('t')
	m_w = Symbol('m_w')
	q1 = Function('q1')
	q2 = Function('q2')
	e = MutableDenseMatrix([[-m_f*s**2*sigma**2*cos(q1(t))/4 - m_f*s**2*sigma**2*cos(q2(t))/4 + 2*m_f*s**2*sigma**2/3 + m_f*s**2*sigma*cos(q1(t))/4 + m_f*s**2*sigma*cos(q2(t))/4 - m_f*s**2*sigma + m_f*s**2/2 + m_w*s**2*sigma**2/12 - m_w*s**2*sigma/6 + m_w*s**2/12, m_f*s**2*sigma**2*cos(q1(t))/8 - m_f*s**2*sigma**2/12 - m_f*s**2*sigma*cos(q1(t))/8, -m_f*s**2*sigma**2*cos(q2(t))/8 + m_f*s**2*sigma**2/12 + m_f*s**2*sigma*cos(q2(t))/8], [m_f*s**2*sigma**2*cos(q1(t))/8 - m_f*s**2*sigma**2/12 - m_f*s**2*sigma*cos(q1(t))/8, m_f*s**2*sigma**2/12, 0], [-m_f*s**2*sigma**2*cos(q2(t))/8 + m_f*s**2*sigma**2/12 + m_f*s**2*sigma*cos(q2(t))/8, 0, m_f*s**2*sigma**2/12]])
	return e
def get_f():
	g = Symbol('g')
	m_f = Symbol('m_f')
	s = Symbol('s')
	sigma = Symbol('sigma')
	t = Symbol('t')
	m_w = Symbol('m_w')
	q0 = Function('q0')
	q1 = Function('q1')
	q2 = Function('q2')
	e = ImmutableDenseMatrix([[-g*m_f*(-s*sigma*(sin(q0(t))*sin(q1(t))/4 + cos(q0(t))*cos(q1(t))/4) - s*(Rational(1, 2) - sigma/2)*cos(q0(t))) - g*m_f*(-s*sigma*(sin(q0(t))*sin(q2(t))/4 - cos(q0(t))*cos(q2(t))/4) - s*(sigma/2 + Rational(-1, 2))*cos(q0(t))) + m_f*s**2*sigma**2*sin(q1(t))*Derivative(q0(t), t)*Derivative(q1(t), t)/4 - m_f*s**2*sigma**2*sin(q1(t))*Derivative(q1(t), t)**2/8 + m_f*s**2*sigma**2*sin(q2(t))*Derivative(q0(t), t)*Derivative(q2(t), t)/4 + m_f*s**2*sigma**2*sin(q2(t))*Derivative(q2(t), t)**2/8 - 2*m_f*s**2*sigma**2*Derivative(q0(t), (t, 2))/3 + 2*m_f*s**2*sigma**2*Derivative(q0(t), (t, 2))/3 - m_f*s**2*sigma**2*Derivative(q1(t), (t, 2))/12 + m_f*s**2*sigma**2*Derivative(q1(t), (t, 2))/12 - m_f*s**2*sigma**2*Derivative(q2(t), (t, 2))/12 + m_f*s**2*sigma**2*Derivative(q2(t), (t, 2))/12 - m_f*s**2*sigma*sin(q1(t))*Derivative(q0(t), t)*Derivative(q1(t), t)/4 + m_f*s**2*sigma*sin(q1(t))*Derivative(q1(t), t)**2/8 - m_f*s**2*sigma*sin(q2(t))*Derivative(q0(t), t)*Derivative(q2(t), t)/4 - m_f*s**2*sigma*sin(q2(t))*Derivative(q2(t), t)**2/8 - m_w*s**2*sigma**2*Derivative(q0(t), (t, 2))/12 + m_w*s**2*sigma**2*Derivative(q0(t), (t, 2))/12], [g*m_f*s*sigma*(-sin(q0(t))*sin(q1(t))/4 - cos(q0(t))*cos(q1(t))/4) - m_f*s**2*sigma**2*sin(q1(t))*Derivative(q0(t), t)*Derivative(q1(t), t)/8 - m_f*s**2*sigma**2*Derivative(q0(t), (t, 2))/12 + m_f*s**2*sigma**2*Derivative(q0(t), (t, 2))/12 - m_f*s**2*sigma**2*Derivative(q1(t), (t, 2))/12 + m_f*s**2*sigma**2*Derivative(q1(t), (t, 2))/12 + m_f*s**2*sigma*sin(q1(t))*Derivative(q0(t), t)*Derivative(q1(t), t)/8 - m_f*s**2*(3*sigma**2*sin(q1(t))*Derivative(q0(t), t)**2 - 3*sigma**2*sin(q1(t))*Derivative(q0(t), t)*Derivative(q1(t), t) - 3*sigma*sin(q1(t))*Derivative(q0(t), t)**2 + 3*sigma*sin(q1(t))*Derivative(q0(t), t)*Derivative(q1(t), t))/24], [g*m_f*s*sigma*(sin(q0(t))*sin(q2(t))/4 - cos(q0(t))*cos(q2(t))/4) + m_f*s**2*sigma**2*sin(q2(t))*Derivative(q0(t), t)*Derivative(q2(t), t)/8 - m_f*s**2*sigma**2*Derivative(q0(t), (t, 2))/12 + m_f*s**2*sigma**2*Derivative(q0(t), (t, 2))/12 - m_f*s**2*sigma**2*Derivative(q2(t), (t, 2))/12 + m_f*s**2*sigma**2*Derivative(q2(t), (t, 2))/12 - m_f*s**2*sigma*sin(q2(t))*Derivative(q0(t), t)*Derivative(q2(t), t)/8 - m_f*s**2*(3*sigma**2*sin(q2(t))*Derivative(q0(t), t)**2 + 3*sigma**2*sin(q2(t))*Derivative(q0(t), t)*Derivative(q2(t), t) - 3*sigma*sin(q2(t))*Derivative(q0(t), t)**2 - 3*sigma*sin(q2(t))*Derivative(q0(t), t)*Derivative(q2(t), t))/24]])
	return e
def get_T():
	m_f = Symbol('m_f')
	s = Symbol('s')
	sigma = Symbol('sigma')
	t = Symbol('t')
	m_w = Symbol('m_w')
	q1 = Function('q1')
	q0 = Function('q0')
	q2 = Function('q2')
	e = m_f*s**2*(-3*sigma**2*cos(q1(t))*Derivative(q0(t), t)**2 + 3*sigma**2*cos(q1(t))*Derivative(q0(t), t)*Derivative(q1(t), t) + 4*sigma**2*Derivative(q0(t), t)**2 - 2*sigma**2*Derivative(q0(t), t)*Derivative(q1(t), t) + sigma**2*Derivative(q1(t), t)**2 + 3*sigma*cos(q1(t))*Derivative(q0(t), t)**2 - 3*sigma*cos(q1(t))*Derivative(q0(t), t)*Derivative(q1(t), t) - 6*sigma*Derivative(q0(t), t)**2 + 3*Derivative(q0(t), t)**2)/24 + m_f*s**2*(-3*sigma**2*cos(q2(t))*Derivative(q0(t), t)**2 - 3*sigma**2*cos(q2(t))*Derivative(q0(t), t)*Derivative(q2(t), t) + 4*sigma**2*Derivative(q0(t), t)**2 + 2*sigma**2*Derivative(q0(t), t)*Derivative(q2(t), t) + sigma**2*Derivative(q2(t), t)**2 + 3*sigma*cos(q2(t))*Derivative(q0(t), t)**2 + 3*sigma*cos(q2(t))*Derivative(q0(t), t)*Derivative(q2(t), t) - 6*sigma*Derivative(q0(t), t)**2 + 3*Derivative(q0(t), t)**2)/24 + m_w*s**2*sigma**2*Derivative(q0(t), t)**2/24 - m_w*s**2*sigma*Derivative(q0(t), t)**2/12 + m_w*s**2*Derivative(q0(t), t)**2/24
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
	a_0 = Symbol('a_0')
	c = Symbol('c')
	rho = Symbol('rho')
	s = Symbol('s')
	sigma = Symbol('sigma')
	t = Symbol('t')
	alpha_1 = Symbol('alpha_1')
	w_g = Symbol('w_g')
	alpha_2 = Symbol('alpha_2')
	q0 = Function('q0')
	q1 = Function('q1')
	q2 = Function('q2')
	e = ImmutableDenseMatrix([[-V*a_0*c*rho*s**3*(1 - sigma)**3*Derivative(q0(t), t)/24 + s**3*sigma**3*(-V*a_0*c*rho*Derivative(q0(t), t)/6 + V*a_0*c*rho*Derivative(q1(t), t)/6)/8 + s**3*sigma**3*(-V*a_0*c*rho*Derivative(q0(t), t)/6 - V*a_0*c*rho*Derivative(q2(t), t)/6)/8 + s**2*sigma**2*(V**2*a_0*alpha_1*c*rho/4 + V*a_0*c*rho*s*sigma*cos(q1(t))*Derivative(q0(t), t)/4 - V*a_0*c*rho*s*sigma*cos(q1(t))*Derivative(q1(t), t)/8 - V*a_0*c*rho*s*cos(q1(t))*Derivative(q0(t), t)/4 + V*a_0*c*rho*s*cos(q1(t))*Derivative(q1(t), t)/8 + V*a_0*c*rho*w_g/4)/4 + s**2*sigma**2*(-V**2*a_0*alpha_2*c*rho/4 + V*a_0*c*rho*s*sigma*cos(q2(t))*Derivative(q0(t), t)/4 + V*a_0*c*rho*s*sigma*cos(q2(t))*Derivative(q2(t), t)/8 - V*a_0*c*rho*s*cos(q2(t))*Derivative(q0(t), t)/4 - V*a_0*c*rho*s*cos(q2(t))*Derivative(q2(t), t)/8 - V*a_0*c*rho*w_g/4)/4 + s*sigma*(-V**2*a_0*alpha_1*c*rho*s*sigma*cos(q1(t))/4 + V**2*a_0*alpha_1*c*rho*s*cos(q1(t))/4 - V*a_0*c*rho*s**2*sigma**2*cos(q1(t))**2*Derivative(q0(t), t)/8 + V*a_0*c*rho*s**2*sigma*cos(q1(t))**2*Derivative(q0(t), t)/4 - V*a_0*c*rho*s**2*cos(q1(t))**2*Derivative(q0(t), t)/8 - V*a_0*c*rho*s*sigma*w_g*cos(q1(t))/4 + V*a_0*c*rho*s*w_g*cos(q1(t))/4)/2 + s*sigma*(V**2*a_0*alpha_2*c*rho*s*sigma*cos(q2(t))/4 - V**2*a_0*alpha_2*c*rho*s*cos(q2(t))/4 - V*a_0*c*rho*s**2*sigma**2*cos(q2(t))**2*Derivative(q0(t), t)/8 + V*a_0*c*rho*s**2*sigma*cos(q2(t))**2*Derivative(q0(t), t)/4 - V*a_0*c*rho*s**2*cos(q2(t))**2*Derivative(q0(t), t)/8 + V*a_0*c*rho*s*sigma*w_g*cos(q2(t))/4 - V*a_0*c*rho*s*w_g*cos(q2(t))/4)/2], [s**3*sigma**3*(V*a_0*c*rho*Derivative(q0(t), t)/6 - V*a_0*c*rho*Derivative(q1(t), t)/6)/8 + s**2*sigma**2*(-V**2*a_0*alpha_1*c*rho/4 - V*a_0*c*rho*s*sigma*cos(q1(t))*Derivative(q0(t), t)/8 + V*a_0*c*rho*s*cos(q1(t))*Derivative(q0(t), t)/8 - V*a_0*c*rho*w_g/4)/4], [s**3*sigma**3*(-V*a_0*c*rho*Derivative(q0(t), t)/6 - V*a_0*c*rho*Derivative(q2(t), t)/6)/8 + s**2*sigma**2*(-V**2*a_0*alpha_2*c*rho/4 + V*a_0*c*rho*s*sigma*cos(q2(t))*Derivative(q0(t), t)/8 - V*a_0*c*rho*s*cos(q2(t))*Derivative(q0(t), t)/8 - V*a_0*c*rho*w_g/4)/4]])
	return e
