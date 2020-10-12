from sympy import *
def get_M():
	I_xxw = Symbol('I_xxw')
	eta_0 = Symbol('eta_0')
	m_w = Symbol('m_w')
	eta_1 = Symbol('eta_1')
	e = MutableDenseMatrix([[I_xxw + eta_0**2*m_w + eta_1**2*m_w]])
	return e
def get_f():
	g = Symbol('g')
	m_w = Symbol('m_w')
	eta_0 = Symbol('eta_0')
	t = Symbol('t')
	eta_1 = Symbol('eta_1')
	q0 = Function('q0')
	e = ImmutableDenseMatrix([[-g*m_w*(-eta_0*cos(q0(t)) + eta_1*sin(q0(t)))]])
	return e
def get_T():
	I_xxw = Symbol('I_xxw')
	t = Symbol('t')
	eta_0 = Symbol('eta_0')
	m_w = Symbol('m_w')
	eta_1 = Symbol('eta_1')
	q0 = Function('q0')
	e = I_xxw*Derivative(q0(t), t)**2/2 + eta_0**2*m_w*Derivative(q0(t), t)**2/2 + eta_1**2*m_w*Derivative(q0(t), t)**2/2
	return e
def get_U():
	g = Symbol('g')
	m_w = Symbol('m_w')
	eta_0 = Symbol('eta_0')
	t = Symbol('t')
	eta_1 = Symbol('eta_1')
	q0 = Function('q0')
	e = g*m_w*(eta_0*sin(q0(t)) + eta_1*cos(q0(t)))
	return e
def get_Q():
	V = Symbol('V')
	a_1 = Symbol('a_1')
	c = Symbol('c')
	rho = Symbol('rho')
	s = Symbol('s')
	sigma = Symbol('sigma')
	alpha_1 = Symbol('alpha_1')
	w_g = Symbol('w_g')
	t = Symbol('t')
	C_Dmax = Symbol('C_Dmax')
	alpha_2 = Symbol('alpha_2')
	a_0 = Symbol('a_0')
	q0 = Function('q0')
	e = ImmutableDenseMatrix([[V**2*a_1*c*rho*s*sigma*(-9*s*sigma/20 + s/2)*(alpha_1 + w_g/V - (-9*s*sigma/20 + s/2)*Derivative(q0(t), t)/V)/20 + V**2*a_1*c*rho*s*sigma*(-7*s*sigma/20 + s/2)*(alpha_1 + w_g/V - (-7*s*sigma/20 + s/2)*Derivative(q0(t), t)/V)/20 + V**2*a_1*c*rho*s*sigma*(-s*sigma/4 + s/2)*(alpha_1 + w_g/V - (-s*sigma/4 + s/2)*Derivative(q0(t), t)/V)/20 + V**2*a_1*c*rho*s*sigma*(-3*s*sigma/20 + s/2)*(alpha_1 + w_g/V - (-3*s*sigma/20 + s/2)*Derivative(q0(t), t)/V)/20 + V**2*a_1*c*rho*s*sigma*(-s*sigma/20 + s/2)*(alpha_1 + w_g/V - (-s*sigma/20 + s/2)*Derivative(q0(t), t)/V)/20 + V**2*c*rho*s*sigma*(s*sigma/20 - s/2)*(C_Dmax*(1 - cos(2*alpha_2 + 2*w_g/V - 2*(s*sigma/20 - s/2)*Derivative(q0(t), t)/V))*sin(alpha_2 - (s*sigma/20 - s/2)*Derivative(q0(t), t)/V)/2 + Float('0.10000000000000001', precision=53)*a_1*(alpha_2 + w_g/V - (s*sigma/20 - s/2)*Derivative(q0(t), t)/V)*cos(alpha_2 - (s*sigma/20 - s/2)*Derivative(q0(t), t)/V))/20 + V**2*c*rho*s*sigma*(3*s*sigma/20 - s/2)*(C_Dmax*(1 - cos(2*alpha_2 + 2*w_g/V - 2*(3*s*sigma/20 - s/2)*Derivative(q0(t), t)/V))*sin(alpha_2 - (3*s*sigma/20 - s/2)*Derivative(q0(t), t)/V)/2 + Float('0.30000000000000004', precision=53)*a_1*(alpha_2 + w_g/V - (3*s*sigma/20 - s/2)*Derivative(q0(t), t)/V)*cos(alpha_2 - (3*s*sigma/20 - s/2)*Derivative(q0(t), t)/V))/20 + V**2*c*rho*s*sigma*(s*sigma/4 - s/2)*(C_Dmax*(1 - cos(2*alpha_2 + 2*w_g/V - 2*(s*sigma/4 - s/2)*Derivative(q0(t), t)/V))*sin(alpha_2 - (s*sigma/4 - s/2)*Derivative(q0(t), t)/V)/2 + Float('0.5', precision=53)*a_1*(alpha_2 + w_g/V - (s*sigma/4 - s/2)*Derivative(q0(t), t)/V)*cos(alpha_2 - (s*sigma/4 - s/2)*Derivative(q0(t), t)/V))/20 + V**2*c*rho*s*sigma*(7*s*sigma/20 - s/2)*(C_Dmax*(1 - cos(2*alpha_2 + 2*w_g/V - 2*(7*s*sigma/20 - s/2)*Derivative(q0(t), t)/V))*sin(alpha_2 - (7*s*sigma/20 - s/2)*Derivative(q0(t), t)/V)/2 + Float('0.70000000000000007', precision=53)*a_1*(alpha_2 + w_g/V - (7*s*sigma/20 - s/2)*Derivative(q0(t), t)/V)*cos(alpha_2 - (7*s*sigma/20 - s/2)*Derivative(q0(t), t)/V))/20 + V**2*c*rho*s*sigma*(9*s*sigma/20 - s/2)*(C_Dmax*(1 - cos(2*alpha_2 + 2*w_g/V - 2*(9*s*sigma/20 - s/2)*Derivative(q0(t), t)/V))*sin(alpha_2 - (9*s*sigma/20 - s/2)*Derivative(q0(t), t)/V)/2 + Float('0.90000000000000002', precision=53)*a_1*(alpha_2 + w_g/V - (9*s*sigma/20 - s/2)*Derivative(q0(t), t)/V)*cos(alpha_2 - (9*s*sigma/20 - s/2)*Derivative(q0(t), t)/V))/20 - V*a_0*c*rho*s**3*(1 - sigma)**3*Derivative(q0(t), t)/24]])
	return e
