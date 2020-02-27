import sympy as sym

class HomogenousTransform:

    def __init__(self,T=None):
        self.T = sym.eye(4) if T is None else T

    def R_x(self,angle):
        H = sym.eye(4)
        H[:3,:3]=sym.Matrix([[1,0,0],
                            [0,sym.cos(angle),-sym.sin(angle)],
                            [0,sym.sin(angle),sym.cos(angle)]])
        return HomogenousTransform(self.T*H)

    def R_y(self,angle):
        H = sym.eye(4)
        H[:3,:3]=sym.Matrix([[sym.cos(angle),0,sym.sin(angle)],
                            [0,1,0],
                            [-sym.sin(angle),0,sym.cos(angle)]])
        return HomogenousTransform(self.T*H)

    def R_z(self,angle):
        H = sym.eye(4)
        H[:3,:3]=sym.Matrix([[sym.cos(angle),-sym.sin(angle),0],
                            [sym.sin(angle),sym.cos(angle),0],
                            [0,0,1]])
        return HomogenousTransform(self.T*H)

    def Translate(self,x,y,z):
        H = sym.eye(4)
        H[:3,3] = sym.Matrix([x,y,z])
        return HomogenousTransform(self.T*H)

    def simplify(self):
        return HomogenousTransform(sym.simplify(self.T))

    def Transform_point(self,p):
        p_l = list(p)
        p_l.append(1)
        p_t = self.T*sym.Matrix(p_l)
        return sym.Matrix(p_t[:3])
