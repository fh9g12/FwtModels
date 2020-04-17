from sympy import *
def get_M():
	c = Symbol('c')
	delta_m = Symbol('delta_m')
	rho_t = Symbol('rho_t')
	s_t = Symbol('s_t')
	sigma = Symbol('sigma')
	eta0 = Symbol('eta0')
	e_0 = Symbol('e_0')
	eta1 = Symbol('eta1')
	e = ImmutableDenseMatrix([[c*delta_m*rho_t*s_t**5*sigma**5/eta0**2 - 4*c*delta_m*rho_t*s_t**5*sigma**4/eta0**2 + 6*c*delta_m*rho_t*s_t**5*sigma**3/eta0**2 - 4*c*delta_m*rho_t*s_t**5*sigma**2/eta0**2 + c*delta_m*rho_t*s_t**5*sigma/eta0**2 - c*rho_t*s_t**5*sigma**5/(5*eta0**2) + c*rho_t*s_t**5*sigma**4/eta0**2 - 2*c*rho_t*s_t**5*sigma**3/eta0**2 + 2*c*rho_t*s_t**5*sigma**2/eta0**2 - c*rho_t*s_t**5*sigma/eta0**2 + c*rho_t*s_t**5/(5*eta0**2), c**2*e_0*rho_t*s_t**4*sigma**4/(4*eta0*eta1) - c**2*e_0*rho_t*s_t**4*sigma**3/(eta0*eta1) + 3*c**2*e_0*rho_t*s_t**4*sigma**2/(2*eta0*eta1) - c**2*e_0*rho_t*s_t**4*sigma/(eta0*eta1) + c**2*e_0*rho_t*s_t**4/(4*eta0*eta1) - c**2*rho_t*s_t**4*sigma**4/(16*eta0*eta1) + c**2*rho_t*s_t**4*sigma**3/(4*eta0*eta1) - 3*c**2*rho_t*s_t**4*sigma**2/(8*eta0*eta1) + c**2*rho_t*s_t**4*sigma/(4*eta0*eta1) - c**2*rho_t*s_t**4/(16*eta0*eta1), -c*delta_m*rho_t*s_t**4*sigma**4/(2*eta0) + c*delta_m*rho_t*s_t**4*sigma**3/eta0 - c*delta_m*rho_t*s_t**4*sigma**2/(2*eta0)], [c**2*e_0*rho_t*s_t**4*sigma**4/(4*eta0*eta1) - c**2*e_0*rho_t*s_t**4*sigma**3/(eta0*eta1) + 3*c**2*e_0*rho_t*s_t**4*sigma**2/(2*eta0*eta1) - c**2*e_0*rho_t*s_t**4*sigma/(eta0*eta1) + c**2*e_0*rho_t*s_t**4/(4*eta0*eta1) - c**2*rho_t*s_t**4*sigma**4/(16*eta0*eta1) + c**2*rho_t*s_t**4*sigma**3/(4*eta0*eta1) - 3*c**2*rho_t*s_t**4*sigma**2/(8*eta0*eta1) + c**2*rho_t*s_t**4*sigma/(4*eta0*eta1) - c**2*rho_t*s_t**4/(16*eta0*eta1), -c**3*e_0**2*rho_t*s_t**3*sigma**3/(3*eta1**2) + c**3*e_0**2*rho_t*s_t**3*sigma**2/eta1**2 - c**3*e_0**2*rho_t*s_t**3*sigma/eta1**2 + c**3*e_0**2*rho_t*s_t**3/(3*eta1**2) + c**3*e_0*rho_t*s_t**3*sigma**3/(6*eta1**2) - c**3*e_0*rho_t*s_t**3*sigma**2/(2*eta1**2) + c**3*e_0*rho_t*s_t**3*sigma/(2*eta1**2) - c**3*e_0*rho_t*s_t**3/(6*eta1**2) - 7*c**3*rho_t*s_t**3*sigma**3/(144*eta1**2) + 7*c**3*rho_t*s_t**3*sigma**2/(48*eta1**2) - 7*c**3*rho_t*s_t**3*sigma/(48*eta1**2) + 7*c**3*rho_t*s_t**3/(144*eta1**2), 0], [-c*delta_m*rho_t*s_t**4*sigma**4/(2*eta0) + c*delta_m*rho_t*s_t**4*sigma**3/eta0 - c*delta_m*rho_t*s_t**4*sigma**2/(2*eta0), 0, c*delta_m*rho_t*s_t**3*sigma**3/3]])
	return e
def get_f():
	EI = Symbol('EI')
	s_t = Symbol('s_t')
	sigma = Symbol('sigma')
	t = Symbol('t')
	eta0 = Symbol('eta0')
	c = Symbol('c')
	delta_m = Symbol('delta_m')
	g = Symbol('g')
	rho_t = Symbol('rho_t')
	alpha_r = Symbol('alpha_r')
	GJ = Symbol('GJ')
	eta1 = Symbol('eta1')
	e_0 = Symbol('e_0')
	q0 = Function('q0')
	q1 = Function('q1')
	e = ImmutableDenseMatrix([[-4*EI*s_t*sigma*q0(t)/eta0**2 + 4*EI*s_t*q0(t)/eta0**2 + c*delta_m*g*rho_t*s_t**3*sigma**3*cos(alpha_r)/eta0 - 2*c*delta_m*g*rho_t*s_t**3*sigma**2*cos(alpha_r)/eta0 + c*delta_m*g*rho_t*s_t**3*sigma*cos(alpha_r)/eta0 - c*g*rho_t*s_t**3*sigma**3*cos(alpha_r)/(3*eta0) + c*g*rho_t*s_t**3*sigma**2*cos(alpha_r)/eta0 - c*g*rho_t*s_t**3*sigma*cos(alpha_r)/eta0 + c*g*rho_t*s_t**3*cos(alpha_r)/(3*eta0)], [-GJ*s_t*sigma*q1(t)/eta1**2 + GJ*s_t*q1(t)/eta1**2 + c**2*e_0*g*rho_t*s_t**2*sigma**2*cos(alpha_r)/(2*eta1) - c**2*e_0*g*rho_t*s_t**2*sigma*cos(alpha_r)/eta1 + c**2*e_0*g*rho_t*s_t**2*cos(alpha_r)/(2*eta1) - c**2*g*rho_t*s_t**2*sigma**2*cos(alpha_r)/(8*eta1) + c**2*g*rho_t*s_t**2*sigma*cos(alpha_r)/(4*eta1) - c**2*g*rho_t*s_t**2*cos(alpha_r)/(8*eta1)], [-c*delta_m*g*rho_t*s_t**2*sigma**2*cos(alpha_r)/2]])
	return e
def get_T():
	e = 0
	return e
def get_U():
	e = 0
	return e
def get_Q():
	V = Symbol('V')
	a_0 = Symbol('a_0')
	alpha_r = Symbol('alpha_r')
	c = Symbol('c')
	rho = Symbol('rho')
	s_t = Symbol('s_t')
	sigma = Symbol('sigma')
	eta0 = Symbol('eta0')
	t = Symbol('t')
	eta1 = Symbol('eta1')
	Lambda = Symbol('Lambda')
	M_thetadot = Symbol('M_thetadot')
	e_0 = Symbol('e_0')
	q1 = Function('q1')
	q0 = Function('q0')
	q2 = Function('q2')
	e = ImmutableDenseMatrix([[V**2*a_0*alpha_r*c*rho*s_t**3*(1 - sigma)**3/(6*eta0) + V**2*a_0*c*rho*s_t**4*(1 - sigma)**4*q1(t)/(8*eta0*eta1) + V**2*a_0*c*rho*s_t**3*sigma*(sigma - 1)**2*atan(sin(alpha_r)/((sin(Lambda)**2 + cos(Lambda)**2)*cos(alpha_r)))/(4*eta0*sqrt(1 + sin(alpha_r)**2/((sin(Lambda)**2 + cos(Lambda)**2)**2*cos(alpha_r)**2))) - V*a_0*c*rho*s_t**5*(1 - sigma)**5*Derivative(q0(t), t)/(10*eta0**2) + (17*V*a_0*c*rho*s_t**4*sigma**2*(sigma - 1)**2/(200*eta0*sqrt(1 + sin(alpha_r)**2/((sin(Lambda)**2 + cos(Lambda)**2)**2*cos(alpha_r)**2))) - 17*V*a_0*c*rho*s_t**4*sigma**2*(sigma - 1)**2*sin(alpha_r)*atan(sin(alpha_r)/((sin(Lambda)**2 + cos(Lambda)**2)*cos(alpha_r)))/(200*eta0*sqrt(1 + sin(alpha_r)**2/((sin(Lambda)**2 + cos(Lambda)**2)**2*cos(alpha_r)**2))*(sin(Lambda)**2 + cos(Lambda)**2)*cos(alpha_r)))*Derivative(q2(t), t) + (-V*a_0*c*rho*s_t**5*sigma*(sigma - 1)**4/(4*eta0**2*sqrt(1 + sin(alpha_r)**2/((sin(Lambda)**2 + cos(Lambda)**2)**2*cos(alpha_r)**2))) + V*a_0*c*rho*s_t**5*sigma*(sigma - 1)**4*sin(alpha_r)*atan(sin(alpha_r)/((sin(Lambda)**2 + cos(Lambda)**2)*cos(alpha_r)))/(4*eta0**2*sqrt(1 + sin(alpha_r)**2/((sin(Lambda)**2 + cos(Lambda)**2)**2*cos(alpha_r)**2))*(sin(Lambda)**2 + cos(Lambda)**2)*cos(alpha_r)))*Derivative(q0(t), t) + (V**2*a_0*c*rho*s_t**3*sigma*(sigma - 1)**2*(sin(Lambda)/(sin(Lambda)**2 + cos(Lambda)**2) + sin(Lambda)*sin(alpha_r)**2/((sin(Lambda)**2 + cos(Lambda)**2)**2*cos(alpha_r)**2))/(4*eta0*(1 + sin(alpha_r)**2/((sin(Lambda)**2 + cos(Lambda)**2)**2*cos(alpha_r)**2))**Rational(3, 2)) - V**2*a_0*c*rho*s_t**3*sigma*(sigma - 1)**2*(sin(Lambda)/(sin(Lambda)**2 + cos(Lambda)**2) + sin(Lambda)*sin(alpha_r)**2/((sin(Lambda)**2 + cos(Lambda)**2)**2*cos(alpha_r)**2))*sin(alpha_r)*atan(sin(alpha_r)/((sin(Lambda)**2 + cos(Lambda)**2)*cos(alpha_r)))/(4*eta0*(1 + sin(alpha_r)**2/((sin(Lambda)**2 + cos(Lambda)**2)**2*cos(alpha_r)**2))**Rational(3, 2)*(sin(Lambda)**2 + cos(Lambda)**2)*cos(alpha_r)))*q2(t) + (V**2*a_0*c*rho*s_t**3*sigma*(sigma - 1)**2*(-2*s_t*(1 - sigma)*sin(Lambda)/(eta0*(sin(Lambda)**2 + cos(Lambda)**2)) - 2*s_t*(1 - sigma)*sin(Lambda)*sin(alpha_r)**2/(eta0*(sin(Lambda)**2 + cos(Lambda)**2)**2*cos(alpha_r)**2))/(4*eta0*(1 + sin(alpha_r)**2/((sin(Lambda)**2 + cos(Lambda)**2)**2*cos(alpha_r)**2))**Rational(3, 2)) - V**2*a_0*c*rho*s_t**3*sigma*(sigma - 1)**2*(-2*s_t*(1 - sigma)*sin(Lambda)/(eta0*(sin(Lambda)**2 + cos(Lambda)**2)) - 2*s_t*(1 - sigma)*sin(Lambda)*sin(alpha_r)**2/(eta0*(sin(Lambda)**2 + cos(Lambda)**2)**2*cos(alpha_r)**2))*sin(alpha_r)*atan(sin(alpha_r)/((sin(Lambda)**2 + cos(Lambda)**2)*cos(alpha_r)))/(4*eta0*(1 + sin(alpha_r)**2/((sin(Lambda)**2 + cos(Lambda)**2)**2*cos(alpha_r)**2))**Rational(3, 2)*(sin(Lambda)**2 + cos(Lambda)**2)*cos(alpha_r)))*q0(t) + (V**2*a_0*c*rho*s_t**4*sigma*(1 - sigma)*(sigma - 1)**2/(4*eta0*eta1*sqrt(1 + sin(alpha_r)**2/((sin(Lambda)**2 + cos(Lambda)**2)**2*cos(alpha_r)**2))) - V**2*a_0*c*rho*s_t**4*sigma*(1 - sigma)*(sigma - 1)**2*sin(alpha_r)*atan(sin(alpha_r)/((sin(Lambda)**2 + cos(Lambda)**2)*cos(alpha_r)))/(4*eta0*eta1*sqrt(1 + sin(alpha_r)**2/((sin(Lambda)**2 + cos(Lambda)**2)**2*cos(alpha_r)**2))*(sin(Lambda)**2 + cos(Lambda)**2)*cos(alpha_r)))*q1(t)], [M_thetadot*V*c**3*rho*s_t**3*(1 - sigma)**3*Derivative(q1(t), t)/(24*eta1**2) + V**2*a_0*alpha_r*c**2*e_0*rho*s_t**2*(1 - sigma)**2/(4*eta1) + V**2*a_0*c**2*e_0*rho*s_t**3*(1 - sigma)**3*q1(t)/(6*eta1**2) - V*a_0*c**2*e_0*rho*s_t**4*(1 - sigma)**4*Derivative(q0(t), t)/(8*eta0*eta1)], [-17*V**2*a_0*c*rho*s_t**2*sigma**2*atan(sin(alpha_r)/((sin(Lambda)**2 + cos(Lambda)**2)*cos(alpha_r)))/(200*sqrt(1 + sin(alpha_r)**2/((sin(Lambda)**2 + cos(Lambda)**2)**2*cos(alpha_r)**2))) + (-17*V*a_0*c*rho*s_t**3*sigma**3/(400*sqrt(1 + sin(alpha_r)**2/((sin(Lambda)**2 + cos(Lambda)**2)**2*cos(alpha_r)**2))) + 17*V*a_0*c*rho*s_t**3*sigma**3*sin(alpha_r)*atan(sin(alpha_r)/((sin(Lambda)**2 + cos(Lambda)**2)*cos(alpha_r)))/(400*sqrt(1 + sin(alpha_r)**2/((sin(Lambda)**2 + cos(Lambda)**2)**2*cos(alpha_r)**2))*(sin(Lambda)**2 + cos(Lambda)**2)*cos(alpha_r)))*Derivative(q2(t), t) + (-17*V**2*a_0*c*rho*s_t**2*sigma**2*(sin(Lambda)/(sin(Lambda)**2 + cos(Lambda)**2) + sin(Lambda)*sin(alpha_r)**2/((sin(Lambda)**2 + cos(Lambda)**2)**2*cos(alpha_r)**2))/(200*(1 + sin(alpha_r)**2/((sin(Lambda)**2 + cos(Lambda)**2)**2*cos(alpha_r)**2))**Rational(3, 2)) + 17*V**2*a_0*c*rho*s_t**2*sigma**2*(sin(Lambda)/(sin(Lambda)**2 + cos(Lambda)**2) + sin(Lambda)*sin(alpha_r)**2/((sin(Lambda)**2 + cos(Lambda)**2)**2*cos(alpha_r)**2))*sin(alpha_r)*atan(sin(alpha_r)/((sin(Lambda)**2 + cos(Lambda)**2)*cos(alpha_r)))/(200*(1 + sin(alpha_r)**2/((sin(Lambda)**2 + cos(Lambda)**2)**2*cos(alpha_r)**2))**Rational(3, 2)*(sin(Lambda)**2 + cos(Lambda)**2)*cos(alpha_r)))*q2(t) + (-17*V**2*a_0*c*rho*s_t**2*sigma**2*(-2*s_t*(1 - sigma)*sin(Lambda)/(eta0*(sin(Lambda)**2 + cos(Lambda)**2)) - 2*s_t*(1 - sigma)*sin(Lambda)*sin(alpha_r)**2/(eta0*(sin(Lambda)**2 + cos(Lambda)**2)**2*cos(alpha_r)**2))/(200*(1 + sin(alpha_r)**2/((sin(Lambda)**2 + cos(Lambda)**2)**2*cos(alpha_r)**2))**Rational(3, 2)) + 17*V**2*a_0*c*rho*s_t**2*sigma**2*(-2*s_t*(1 - sigma)*sin(Lambda)/(eta0*(sin(Lambda)**2 + cos(Lambda)**2)) - 2*s_t*(1 - sigma)*sin(Lambda)*sin(alpha_r)**2/(eta0*(sin(Lambda)**2 + cos(Lambda)**2)**2*cos(alpha_r)**2))*sin(alpha_r)*atan(sin(alpha_r)/((sin(Lambda)**2 + cos(Lambda)**2)*cos(alpha_r)))/(200*(1 + sin(alpha_r)**2/((sin(Lambda)**2 + cos(Lambda)**2)**2*cos(alpha_r)**2))**Rational(3, 2)*(sin(Lambda)**2 + cos(Lambda)**2)*cos(alpha_r)))*q0(t) + (17*V*a_0*c*rho*s_t**4*sigma**2*(sigma - 1)**2/(200*eta0*sqrt(1 + sin(alpha_r)**2/((sin(Lambda)**2 + cos(Lambda)**2)**2*cos(alpha_r)**2))) - 17*V*a_0*c*rho*s_t**4*sigma**2*(sigma - 1)**2*sin(alpha_r)*atan(sin(alpha_r)/((sin(Lambda)**2 + cos(Lambda)**2)*cos(alpha_r)))/(200*eta0*sqrt(1 + sin(alpha_r)**2/((sin(Lambda)**2 + cos(Lambda)**2)**2*cos(alpha_r)**2))*(sin(Lambda)**2 + cos(Lambda)**2)*cos(alpha_r)))*Derivative(q0(t), t) + (-17*V**2*a_0*c*rho*s_t**3*sigma**2*(1 - sigma)/(200*eta1*sqrt(1 + sin(alpha_r)**2/((sin(Lambda)**2 + cos(Lambda)**2)**2*cos(alpha_r)**2))) + 17*V**2*a_0*c*rho*s_t**3*sigma**2*(1 - sigma)*sin(alpha_r)*atan(sin(alpha_r)/((sin(Lambda)**2 + cos(Lambda)**2)*cos(alpha_r)))/(200*eta1*sqrt(1 + sin(alpha_r)**2/((sin(Lambda)**2 + cos(Lambda)**2)**2*cos(alpha_r)**2))*(sin(Lambda)**2 + cos(Lambda)**2)*cos(alpha_r)))*q1(t)]])
	return e
