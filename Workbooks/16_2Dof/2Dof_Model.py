from sympy import *
import moyra as ma

def get_M(p):
	t = Symbol('t')
	e = MutableDenseMatrix([[p.m + p.m_w, -p.l_com*p.m*cos(p.q[1])], [-p.l_com*p.m*cos(p.q[1]), p.I_xx + p.l_com**2*p.m]])
	return e
def get_f(p):
	t = Symbol('t')
	e = ImmutableDenseMatrix([[p.g*p.m + p.g*p.m_w + p.k_w*p.q[0] + p.l_com*p.m*sin(p.q[1])*p.qd[1]**2], [-p.g*p.l_com*p.m*cos(p.q[1])]])
	return e
def get_T(p):
	t = Symbol('t')
	e = p.I_xx*p.qd[1]**2/2 + p.l_com**2*p.m*p.qd[1]**2/2 - p.l_com*p.m*cos(p.q[1])*p.qd[0]*p.qd[1] + p.m*p.qd[0]**2/2 + p.m_w*p.qd[0]**2/2
	return e
def get_U(p):
	t = Symbol('t')
	e = p.g*p.m*(-p.l_com*sin(p.q[1]) + p.q[0]) + p.g*p.m_w*p.q[0] + p.k_w*p.q[0]**2/2
	return e
def get_Q(p):
	t = Symbol('t')
	e = MutableDenseMatrix([[p.a[0]*p.V**2*p.c*p.rho*p.s*(p.alpha_1 - (-p.s*p.qd[1]/20 + cos(p.q[1])*p.qd[0])/p.V)*cos(p.q[1])/20 + p.a[1]*p.V**2*p.c*p.rho*p.s*(p.alpha_1 - (-3*p.s*p.qd[1]/20 + cos(p.q[1])*p.qd[0])/p.V)*cos(p.q[1])/20 + p.a[2]*p.V**2*p.c*p.rho*p.s*(p.alpha_1 - (-p.s*p.qd[1]/4 + cos(p.q[1])*p.qd[0])/p.V)*cos(p.q[1])/20 + p.a[3]*p.V**2*p.c*p.rho*p.s*(p.alpha_1 - (-7*p.s*p.qd[1]/20 + cos(p.q[1])*p.qd[0])/p.V)*cos(p.q[1])/20 + p.a[4]*p.V**2*p.c*p.rho*p.s*(p.alpha_1 - (-9*p.s*p.qd[1]/20 + cos(p.q[1])*p.qd[0])/p.V)*cos(p.q[1])/20 + p.a[5]*p.V**2*p.c*p.rho*p.s*(p.alpha_1 - (-11*p.s*p.qd[1]/20 + cos(p.q[1])*p.qd[0])/p.V)*cos(p.q[1])/20 + p.a[6]*p.V**2*p.c*p.rho*p.s*(p.alpha_1 - (-13*p.s*p.qd[1]/20 + cos(p.q[1])*p.qd[0])/p.V)*cos(p.q[1])/20 + p.a[7]*p.V**2*p.c*p.rho*p.s*(p.alpha_1 - (-3*p.s*p.qd[1]/4 + cos(p.q[1])*p.qd[0])/p.V)*cos(p.q[1])/20 + p.a[8]*p.V**2*p.c*p.rho*p.s*(p.alpha_1 - (-17*p.s*p.qd[1]/20 + cos(p.q[1])*p.qd[0])/p.V)*cos(p.q[1])/20 + p.a[9]*p.V**2*p.c*p.rho*p.s*(p.alpha_1 - (-19*p.s*p.qd[1]/20 + cos(p.q[1])*p.qd[0])/p.V)*cos(p.q[1])/20], [-p.a[0]*p.V**2*p.c*p.rho*p.s**2*(p.alpha_1 - (-p.s*p.qd[1]/20 + cos(p.q[1])*p.qd[0])/p.V)/400 - 3*p.a[1]*p.V**2*p.c*p.rho*p.s**2*(p.alpha_1 - (-3*p.s*p.qd[1]/20 + cos(p.q[1])*p.qd[0])/p.V)/400 - p.a[2]*p.V**2*p.c*p.rho*p.s**2*(p.alpha_1 - (-p.s*p.qd[1]/4 + cos(p.q[1])*p.qd[0])/p.V)/80 - 7*p.a[3]*p.V**2*p.c*p.rho*p.s**2*(p.alpha_1 - (-7*p.s*p.qd[1]/20 + cos(p.q[1])*p.qd[0])/p.V)/400 - 9*p.a[4]*p.V**2*p.c*p.rho*p.s**2*(p.alpha_1 - (-9*p.s*p.qd[1]/20 + cos(p.q[1])*p.qd[0])/p.V)/400 - 11*p.a[5]*p.V**2*p.c*p.rho*p.s**2*(p.alpha_1 - (-11*p.s*p.qd[1]/20 + cos(p.q[1])*p.qd[0])/p.V)/400 - 13*p.a[6]*p.V**2*p.c*p.rho*p.s**2*(p.alpha_1 - (-13*p.s*p.qd[1]/20 + cos(p.q[1])*p.qd[0])/p.V)/400 - 3*p.a[7]*p.V**2*p.c*p.rho*p.s**2*(p.alpha_1 - (-3*p.s*p.qd[1]/4 + cos(p.q[1])*p.qd[0])/p.V)/80 - 17*p.a[8]*p.V**2*p.c*p.rho*p.s**2*(p.alpha_1 - (-17*p.s*p.qd[1]/20 + cos(p.q[1])*p.qd[0])/p.V)/400 - 19*p.a[9]*p.V**2*p.c*p.rho*p.s**2*(p.alpha_1 - (-19*p.s*p.qd[1]/20 + cos(p.q[1])*p.qd[0])/p.V)/400]])
	return e
def get_p():
	p = ma.DynamicModelParameters(2)
	p.c = ma.ModelSymbol(value=0.15, string='c')
	p.m = ma.ModelSymbol(value=1, string='m')
	p.m_1 = ma.ModelSymbol(value=1, string='m_1')
	p.m_w = ma.ModelSymbol(value=4, string='m_w')
	p.s = ma.ModelSymbol(value=1, string='s')
	p.f_0 = ma.ModelSymbol(value=2, string='f_0')
	p.k_w = ma.ModelSymbol(value=2, string='k_w')
	p.c_w = ma.ModelSymbol(value=2, string='c_w')
	p.d_w = ma.ModelSymbol(value=2, string='d_w')
	p.d_a = ma.ModelSymbol(value=2, string='d_a')
	p.k_fwt = ma.ModelSymbol(value=0, string='k_fwt')
	p.I_xx = ma.ModelSymbol(value=1, string='I_xx')
	p.l_com = ma.ModelSymbol(value=0.5, string='l_com')
	p.l_m = ma.ModelSymbol(value=0.5, string='l_m')
	p.Lambda = ma.ModelSymbol(value=0.17453292519943295, string='Lambda')
	p.rho = ma.ModelSymbol(value=1.225, string='rho')
	p.V = ma.ModelSymbol(value=10, string='V')
	p.g = ma.ModelSymbol(value=9.81, string='g')
	p.alpha_r = ma.ModelSymbol(value=0.05235987755982989, string='alpha_r')
	p.M_thetadot = ma.ModelSymbol(value=1.2, string='M_thetadot')
	p.alpha_1 = ma.ModelSymbol(value=0, string='alpha_1')
	p.alphadot_1 = ma.ModelSymbol(value=0, string='alphadot_1')
	p.clip_factor = ma.ModelSymbol(value=100, string='mu')
	p.c_d_max = ma.ModelSymbol(value=1, string='C_Dmax')
	p.w_g = ma.ModelSymbol(value=0, string='w_g')
	p.y_0 = Symbol('y_0')
	p.y_i = Symbol('y_i')
	p.x_0 = Symbol('x_0')
	p.fp = ma.ModelMatrix(value=[0, 0, 0, 0], string='qtilde', length=4)
	p.a = ma.ModelMatrix(value=[6.283185307179586, 6.283185307179586, 6.283185307179586, 6.283185307179586, 6.283185307179586, 6.283185307179586, 6.283185307179586, 6.283185307179586, 6.283185307179586, 6.283185307179586], string='a', length=10)
	p.a_0 = ma.ModelSymbol(value=6.283185307179586, string='a_0')
	return p
