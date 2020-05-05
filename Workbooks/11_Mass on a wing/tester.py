from sympy import *
def get_M():
	c = Symbol('c')
	rho = Symbol('rho')
	s = Symbol('s')
	eta0 = Symbol('eta0')
	m_1 = Symbol('m_1')
	eta1 = Symbol('eta1')
	eta2 = Symbol('eta2')
	e = MutableDenseMatrix([[c*rho*s**5*eta0**2/5 + m_1*s**4*eta0**2, c*rho*s**6*eta0*eta1/6 + m_1*s**5*eta0*eta1, c*rho*s**7*eta0*eta2/7 + m_1*s**6*eta0*eta2], [c*rho*s**6*eta0*eta1/6 + m_1*s**5*eta0*eta1, c*rho*s**7*eta1**2/7 + m_1*s**6*eta1**2, c*rho*s**8*eta1*eta2/8 + m_1*s**7*eta1*eta2], [c*rho*s**7*eta0*eta2/7 + m_1*s**6*eta0*eta2, c*rho*s**8*eta1*eta2/8 + m_1*s**7*eta1*eta2, c*rho*s**9*eta2**2/9 + m_1*s**8*eta2**2]])
	return e
def get_f():
	EI = Symbol('EI')
	s = Symbol('s')
	eta0 = Symbol('eta0')
	eta2 = Symbol('eta2')
	t = Symbol('t')
	eta1 = Symbol('eta1')
	c = Symbol('c')
	g = Symbol('g')
	rho = Symbol('rho')
	m_1 = Symbol('m_1')
	q2 = Function('q2')
	q1 = Function('q1')
	q0 = Function('q0')
	e = ImmutableDenseMatrix([[8*EI*s**3*eta0*eta2*q2(t) + 6*EI*s**2*eta0*eta1*q1(t) + 4*EI*s*eta0**2*q0(t) + c*g*rho*s**3*eta0/3 + g*m_1*s**2*eta0], [18*EI*s**4*eta1*eta2*q2(t) + 12*EI*s**3*eta1**2*q1(t) + 6*EI*s**2*eta0*eta1*q0(t) + c*g*rho*s**4*eta1/4 + g*m_1*s**3*eta1], [144*EI*s**5*eta2**2*q2(t)/5 + 18*EI*s**4*eta1*eta2*q1(t) + 8*EI*s**3*eta0*eta2*q0(t) + c*g*rho*s**5*eta2/5 + g*m_1*s**4*eta2]])
	return e
def get_T():
	c = Symbol('c')
	rho = Symbol('rho')
	s = Symbol('s')
	eta2 = Symbol('eta2')
	t = Symbol('t')
	eta1 = Symbol('eta1')
	eta0 = Symbol('eta0')
	m_1 = Symbol('m_1')
	q2 = Function('q2')
	q1 = Function('q1')
	q0 = Function('q0')
	e = c*rho*s**9*eta2**2*Derivative(q2(t), t)**2/18 + c*rho*s**8*eta1*eta2*Derivative(q1(t), t)*Derivative(q2(t), t)/8 + c*rho*s**7*eta0*eta2*Derivative(q0(t), t)*Derivative(q2(t), t)/7 + c*rho*s**7*eta1**2*Derivative(q1(t), t)**2/14 + c*rho*s**6*eta0*eta1*Derivative(q0(t), t)*Derivative(q1(t), t)/6 + c*rho*s**5*eta0**2*Derivative(q0(t), t)**2/10 + m_1*s**8*eta2**2*Derivative(q2(t), t)**2/2 + m_1*s**7*eta1*eta2*Derivative(q1(t), t)*Derivative(q2(t), t) + m_1*s**6*eta0*eta2*Derivative(q0(t), t)*Derivative(q2(t), t) + m_1*s**6*eta1**2*Derivative(q1(t), t)**2/2 + m_1*s**5*eta0*eta1*Derivative(q0(t), t)*Derivative(q1(t), t) + m_1*s**4*eta0**2*Derivative(q0(t), t)**2/2
	return e
def get_U():
	EI = Symbol('EI')
	s = Symbol('s')
	eta2 = Symbol('eta2')
	t = Symbol('t')
	eta1 = Symbol('eta1')
	eta0 = Symbol('eta0')
	c = Symbol('c')
	g = Symbol('g')
	rho = Symbol('rho')
	m_1 = Symbol('m_1')
	q2 = Function('q2')
	q1 = Function('q1')
	q0 = Function('q0')
	e = 72*EI*s**5*eta2**2*q2(t)**2/5 + 18*EI*s**4*eta1*eta2*q1(t)*q2(t) + 6*EI*s**2*eta0*eta1*q0(t)*q1(t) + 2*EI*s*eta0**2*q0(t)**2 + c*g*rho*s**5*eta2*q2(t)/5 + c*g*rho*s**4*eta1*q1(t)/4 + c*g*rho*s**3*eta0*q0(t)/3 + g*m_1*(s**4*eta2*q2(t) + s**3*eta1*q1(t) + s**2*eta0*q0(t)) + s**3*(8*EI*eta0*eta2*q0(t)*q2(t) + 6*EI*eta1**2*q1(t)**2)
	return e
def get_Q():
	sigma = Symbol('sigma')
	t = Symbol('t')
	q0 = Function('q0')
	q1 = Function('q1')
	q2 = Function('q2')
	e = MutableDenseMatrix([[-sigma*Derivative(q0(t), t)], [-sigma*Derivative(q1(t), t)], [-sigma*Derivative(q2(t), t)]])
	return e
