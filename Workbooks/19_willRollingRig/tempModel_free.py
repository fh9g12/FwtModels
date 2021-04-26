from sympy import *
def get_M():
	I_xxf = Symbol('I_xxf')
	I_xxw = Symbol('I_xxw')
	l_f = Symbol('l_f')
	m_f = Symbol('m_f')
	s = Symbol('s')
	sigma = Symbol('sigma')
	m_w = Symbol('m_w')
	y_w = Symbol('y_w')
	z_w = Symbol('z_w')
	e = ImmutableDenseMatrix([[2*I_xxf + I_xxw + 2*l_f**2*m_f - 2*l_f*m_f*s*sigma + 2*l_f*m_f*s + m_f*s**2*sigma**2/2 - m_f*s**2*sigma + m_f*s**2/2 + m_w*y_w**2 + m_w*z_w**2, -I_xxf - l_f**2*m_f + l_f*m_f*s*sigma/2 - l_f*m_f*s/2, I_xxf + l_f**2*m_f - l_f*m_f*s*sigma/2 + l_f*m_f*s/2], [-I_xxf - l_f**2*m_f + l_f*m_f*s*sigma/2 - l_f*m_f*s/2, I_xxf + l_f**2*m_f, 0], [I_xxf + l_f**2*m_f - l_f*m_f*s*sigma/2 + l_f*m_f*s/2, 0, I_xxf + l_f**2*m_f]])
	return e
def get_f():
	g = Symbol('g')
	m_w = Symbol('m_w')
	y_w = Symbol('y_w')
	z_w = Symbol('z_w')
	t = Symbol('t')
	l_f = Symbol('l_f')
	m_f = Symbol('m_f')
	q0 = Function('q0')
	e = ImmutableDenseMatrix([[-g*m_w*y_w + g*m_w*z_w*q0(t)], [g*l_f*m_f], [g*l_f*m_f]])
	return e
def get_T():
	e = 0
	return e
def get_U():
	e = 0
	return e
def get_Q():
	V = Symbol('V')
	a0 = Symbol('a0')
	alpha_c = Symbol('alpha_c')
	c = Symbol('c')
	rho = Symbol('rho')
	s = Symbol('s')
	sigma = Symbol('sigma')
	t = Symbol('t')
	Lambda = Symbol('Lambda')
	a3 = Symbol('a3')
	a1 = Symbol('a1')
	a2 = Symbol('a2')
	q2 = Function('q2')
	q1 = Function('q1')
	q0 = Function('q0')
	e = MutableDenseMatrix([[-V**2*a0*alpha_c*c*rho*s*sigma*(s*sigma/4 - s/2)/4 + V**2*a0*c*rho*s*sigma*(s*sigma/4 - s/2)*q2(t)*sin(Lambda)/4 - V**2*a3*alpha_c*c*rho*s*sigma*(-s*sigma/4 + s/2)/4 + V**2*a3*c*rho*s*sigma*(-s*sigma/4 + s/2)*q1(t)*sin(Lambda)/4 + V*a0*c*rho*s**2*sigma**2*(s*sigma/4 - s/2)*Derivative(q2(t), t)/16 + V*a3*c*rho*s**2*sigma**2*(-s*sigma/4 + s/2)*Derivative(q1(t), t)/16 + (-V*a0*c*rho*s*sigma*(s*sigma/4 - s/2)**2/4 - V*a1*c*rho*s**3*(1 - sigma)**3/64 - V*a2*c*rho*s**3*(1 - sigma)**3/64 - V*a3*c*rho*s*sigma*(-s*sigma/4 + s/2)**2/4)*Derivative(q0(t), t)], [V**2*a3*alpha_c*c*rho*s**2*sigma**2/16 - V**2*a3*c*rho*s**2*sigma**2*q1(t)*sin(Lambda)/16 - V*a3*c*rho*s**3*sigma**3*Derivative(q1(t), t)/64 + V*a3*c*rho*s**2*sigma**2*(-s*sigma/4 + s/2)*Derivative(q0(t), t)/16], [V**2*a0*alpha_c*c*rho*s**2*sigma**2/16 - V**2*a0*c*rho*s**2*sigma**2*q2(t)*sin(Lambda)/16 - V*a0*c*rho*s**3*sigma**3*Derivative(q2(t), t)/64 + V*a0*c*rho*s**2*sigma**2*(s*sigma/4 - s/2)*Derivative(q0(t), t)/16]])
	return e
