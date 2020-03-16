import sympy as sym
from sympy.abc import t

def Vee(E):
    val = sym.Matrix([[0]]*6)
    val[:3,0] = E[:3,3]
    val[3,0] = E[2,1]
    val[4,0] = E[0,2]
    val[5,0] = E[1,0]
    return val

def Wedge(V):
    val = sym.Matrix([[0]*4]*4)
    val[:3,3] = V[:3,0]   
    val[0,1] = -V[-1]
    val[1,0] = V[-1]
    val[2,0] = -V[-2]
    val[0,2] = V[-2]
    val[1,2] = -V[-3]
    val[2,1] = V[-3]
    return val

def Wedge3(V):
    val = sym.Matrix([[0]*3]*3)
    val[0,1] = -V[2]
    val[1,0] = V[2]
    val[2,0] = -V[1]
    val[0,2] = V[1]
    val[1,2] = -V[0]
    val[2,1] = V[0]
    return val


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

    def ManipJacobian(self,q):
        inv = self.Inverse().E
        J = sym.zeros(6,len(q))
        for i,qi in enumerate(q):
            J[:,i] = Vee(self.E.diff(qi)*inv)
        return sym.simplify(J)

    def Adjoint(self):
        E = sym.zeros(6,6)
        E[:3,:3] = self.R
        E[3:,3:] = self.R
        E[:3,3:] = Wedge3(self.t)*self.R.T
        return E

    def InvAdjoint(self):
        E = sym.zeros(6,6)
        E[:3,:3] = self.R.T
        E[3:,3:] = self.R.T
        E[:3,3:] = -self.R.T*Wedge3(self.t)
        return E

    def diff(self,d):
        return HomogenousTransform(self.E.diff(d))

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

    def R_x_small(self,angle):
        H = sym.eye(4)
        H[:3,:3]=sym.Matrix([[1,0,0],
                            [0,1,-angle],
                            [0,angle,1]])
        return HomogenousTransform(self.E*H)

    def R_y_small(self,angle):
        H = sym.eye(4)
        H[:3,:3]=sym.Matrix([[1,0,angle],
                            [0,1,0],
                            [-angle,0,1]])
        return HomogenousTransform(self.E*H)

    def R_z_small(self,angle):
        H = sym.eye(4)
        H[:3,:3]=sym.Matrix([[1,-angle,0],
                            [angle,1,0],
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
