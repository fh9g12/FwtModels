from sympy import *
def get_aoa():
	Lambda = Symbol('Lambda')
	t = Symbol('t')
	alpha_r = Symbol('alpha_r')
	w_g = Symbol('w_g')
	V = Symbol('V')
	q1 = Function('q1')
	e = atan((sin(Lambda)*sin(q1(t))*cos(alpha_r + w_g/V) + sin(alpha_r + w_g/V)*cos(q1(t)))/((sin(Lambda)**2*cos(q1(t)) + cos(Lambda)**2)*cos(alpha_r + w_g/V) - sin(Lambda)*sin(alpha_r + w_g/V)*sin(q1(t))))
	return e
