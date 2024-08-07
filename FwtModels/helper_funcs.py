import sympy as sym

def ShapeFunctions_BN_TM(n,m,q,y_s,x,x_f,alpha_r,factor = 1):
    # check q is the length of n+m
    if n+m != len(q):
        raise ValueError('the sum of n+m must be the same as a length of q')

    # make factor a list the size of n+m
    if isinstance(factor,int) | isinstance(factor,float):
        factor = [factor]*(n+m)
    z = sym.Integer(0)
    tau = alpha_r

    for i in range(0,n):
        z = z + q[i]*y_s**(2+i)*factor[i]
    for i in range(0,m):
        qi = i+n
        tau = tau + q[qi]*y_s**(i+1)*factor[n+i]
    
    z -= tau*(x-x_f)

    return z, tau

def GetVh(alpha,beta,Lambda,theta,mu,simplify=True):
    ## a/c to wind transform
    Ac_V = sym.rot_axis3(beta)*sym.rot_axis2(alpha)   # transform from a/c to the velocity frame 
    ## a/c to Hinge transform
    H_Ac = sym.rot_axis3(-Lambda)*\
                    sym.rot_axis1(-theta)*  \
                    sym.rot_axis3(Lambda)*sym.rot_axis1(-mu)       # transform from a/c to the hinge frame 
    if simplify:
        H_Ac = sym.trigsimp(H_Ac)

    H_V = H_Ac*Ac_V
    if simplify:
        H_V = sym.trigsimp(H_V) # transform from velocity to hinge reference frame

    # Velocity vector in velocity frame is of the form [v 0 0]
    Vv = sym.Matrix([1,0,0])
    # Transform into the hinge reference frame
    return (H_V * Vv)

def GetAoA(alpha,beta,Lambda,theta,mu,simplify=True):
    Vh = GetVh(alpha,beta,Lambda,theta,mu,simplify)
    return sym.atan(Vh[2]/Vh[0])