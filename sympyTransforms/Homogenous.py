import sympy as sym
from sympy.abc import t

class HomogenousTransform:

    def __init__(self,T=None):
        self.E = sym.eye(4) if T is None else T
        self.R = self.E[:3,:3]
        self.t = self.E[:3,3]

    def Inverse(self):
        E = sym.eye(4)
        E[:3,:3] = self.R.T
        E[:3,3] = -self.R.T*self.t
        return HomogenousTransform(E)

    def BodyVelocity(self):
        V = sym.ones(6,1)
        V[:3,0] = self.R.T*self.t.diff(t)

        # Angular velocities skew symetric matrix
        S = self.R.T*self.R.diff(t)
        V[3,0] = S[2,1]
        V[4,0] = S[0,2]
        V[5,0] = S[1,0]
        return V

    def R_x(self,angle):
        H = sym.eye(4)
        H[:3,:3]=sym.Matrix([[1,0,0],
                            [0,sym.cos(angle),-sym.sin(angle)],
                            [0,sym.sin(angle),sym.cos(angle)]])
        return HomogenousTransform(self.E*H)

    def R_y(self,angle):
        H = sym.eye(4)
        H[:3,:3]=sym.Matrix([[sym.cos(angle),0,sym.sin(angle)],
                            [0,1,0],
                            [-sym.sin(angle),0,sym.cos(angle)]])
        return HomogenousTransform(self.E*H)

    def R_z(self,angle):
        H = sym.eye(4)
        H[:3,:3]=sym.Matrix([[sym.cos(angle),-sym.sin(angle),0],
                            [sym.sin(angle),sym.cos(angle),0],
                            [0,0,1]])
        return HomogenousTransform(self.E*H)

    def Translate(self,x,y,z):
        H = sym.eye(4)
        H[:3,3] = sym.Matrix([x,y,z])
        return HomogenousTransform(self.E*H)

    def simplify(self):
        return HomogenousTransform(sym.simplify(self.E))

    def Transform_point(self,p):
        p_l = list(p)
        p_l.append(1)
        p_t = self.E*sym.Matrix(p_l)
        return sym.Matrix(p_t[:3])
