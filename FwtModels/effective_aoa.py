import sympy as sym

def GeomExact(alpha,beta,Lambda,theta,V):
    ## a/c to wind transform
    Ac_V = sym.rot_axis3(beta)*sym.rot_axis2(alpha)

    ## a/c to Hinge transform
    H_Ac = sym.trigsimp(sym.rot_axis3(-Lambda)*\
                        sym.rot_axis1(theta)*  \
                        sym.rot_axis3(Lambda))        # transform from a/c to the hinge frame 


    H_V = sym.trigsimp(H_Ac*Ac_V)                      # transform from velocity to hinge reference frame

    # Transform into the hinge reference frame
    Vh = H_V * V
    return Vh

