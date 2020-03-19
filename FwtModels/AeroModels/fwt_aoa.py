import sympy as sym

def GetVh(alpha,beta,Lambda,theta):
    ## a/c to wind transform
    Ac_V = sym.rot_axis3(beta)*sym.rot_axis2(alpha)   # transform from a/c to the velocity frame 
    ## a/c to Hinge transform
    H_Ac = sym.trigsimp(sym.rot_axis3(-Lambda)*\
                    sym.rot_axis1(theta)*  \
                    sym.rot_axis3(Lambda))        # transform from a/c to the hinge frame 


    H_V = sym.trigsimp(H_Ac*Ac_V) # transform from velocity to hinge reference frame

    # Velocity vector in velocity frame is of the form [v 0 0]
    Vv = sym.Matrix([1,0,0])
    # Transform into the hinge reference frame
    return sym.simplify(H_V * Vv)

def GetAoA(alpha,beta,Lambda,theta):
    Vh = GetVh(alpha,beta,Lambda,theta)
    return sym.atan(Vh[2]/Vh[0])