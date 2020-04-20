import sympy as sym
from sympy.abc import t
import numpy as np
import sys, os

sys.path.insert(1, os.path.join(sys.path[0], '../..'))

import ModelFramework as mf
import ModelFramework.Elements as ele
import ModelFramework.ExternalForces as ef
import FwtModels.RectWing as rw
import sympy.physics.mechanics as me

def GenRectWingModel(b_modes,t_modes,fwt_free,iwt,iwb,fwt_frot):
    p = rw.base_params(b_modes+t_modes+1)

    #get shape functions for main wing
    z_0,tau_0 = mf.ShapeFunctions_BN_TM(b_modes,t_modes,p.q[:-1],p.y_0,p.x_0,p.x_f0,0,factor=p.eta)

    #define wrefernce frames
    wing_root_frame = mf.HomogenousTransform().R_y(p.alpha_r)
    wing_frame = wing_root_frame.Translate(p.x_0,p.y_0,z_0)
    wing_flexural_frame = wing_frame.msubs({p.x_0:p.x_f0})

    fwt_root_frame = wing_frame.msubs({p.y_0:p.s_0,p.x_0:p.x_f0}).Translate(-p.x_f0,0,0)
    if fwt_free:
        fwt_root_frame = fwt_root_frame.R_x(-p.q[-1])
    fwt_flexural_frame = fwt_root_frame.Translate(p.x_f1,p.y_1,0)
    fwt_com_frame = fwt_root_frame.Translate(p.c/2,p.s_1/2,0)
    
    #Create Elemnts
    M_wing = ele.MassMatrix(p.rho_t)
    M_fwt = ele.MassMatrix(p.m_1,I_xx = p.I_xx_1)

    inner_wing_ele = ele.FlexiElement(wing_root_frame,M_wing,p.x_0,p.y_0,z_0,p.c,p.s_0,p.x_f0,p.EI,p.GJ,gravityPot=True)
    fwt_ele = ele.RigidElement(fwt_com_frame,M_fwt,True)

    # Create AeroForces
    wing_AeroForces = ef.AeroForce_1.PerUnitSpan(p,wing_flexural_frame,p.a_0,
                                alphadot = tau_0 if isinstance(tau_0,int) else tau_0.diff(t),
                                M_thetadot = p.M_thetadot,
                                e = p.e_0,
                                rootAlpha = p.alpha_r,
                                deltaAlpha = tau_0,
                                alpha_zero = 0).integrate((p.y_0,0,p.s_0))
  
    alpha_fwt =  0
    alphadot_fwt = 0

    if fwt_free:
        alpha_fwt +=  p.alpha_1
        alphadot_fwt += p.alphadot_1
    #if iwt:
    #    tau_s0 = tau_0.subs(p.y_0,p.s_0)
    #    alpha_fwt += tau_s0
    #    alphadot_fwt += tau_s0.diff(t)

    if fwt_frot:
        fwt_AeroForces_perUnit = ef.AeroForce_2.PerUnitSpan(p,fwt_flexural_frame,p.a_1,
                            alphadot = p.alphadot_1,
                            M_thetadot = p.M_thetadot,
                            e = p.e_1,
                            rootAlpha = p.alpha_1,
                            alpha_zero = 0)
    else:
        fwt_AeroForces_perUnit = ef.AeroForce_1.PerUnitSpan(p,fwt_flexural_frame,p.a_1,
                            alphadot = p.alphadot_1,
                            M_thetadot = p.M_thetadot,
                            e = p.e_1,
                            rootAlpha = p.alpha_1,
                            deltaAlpha = 0, 
                            alpha_zero = 0)

    forces = []
    segments = 5
    for i in range(segments):
        seg_width = p.s_1/segments
        yi = seg_width/2 + i*seg_width
        forces.append(fwt_AeroForces_perUnit.subs({p.y_1:yi})*seg_width)
    Q = sym.Matrix([0]*p.qs)
    for f in forces:
        Q += f.Q()
    fwt_AeroForces = ef.ExternalForce(Q)

    # Setup AoA of FWT
    fwt_aoa = mf.GetAoA(p.alpha_r,0,p.Lambda,0 if not fwt_free else p.q[-1])

    if iwb:
        wing_bend = sym.atan(z_0.diff(p.y_0).subs({p.x_0:p.x_f0,p.y_0:p.s_0}))
        fwt_aoa = me.msubs(fwt_aoa,{p.q[-1]:p.q[-1]-wing_bend})
    if iwt:
        tau_s0 = tau_0.subs(p.y_0,p.s_0)
        fwt_aoa = me.msubs(fwt_aoa,{p.alpha_r:p.alpha_r+tau_s0})

    ## Sub in Aero Forces
    fwt_AeroForces = fwt_AeroForces.subs({p.alpha_1:fwt_aoa,p.alphadot_1:fwt_aoa.diff(t)})

    #Create Composite force
    CompositeForce = ef.CompositeForce([wing_AeroForces,fwt_AeroForces])
    # Create the SYmbolic Model
    sm = mf.SymbolicModel.FromElementsAndForces(p,[inner_wing_ele,fwt_ele],CompositeForce)

    return sm,p


def Gen2DofModel(fwt_free,fwt_frot,rot_AoA=True):
    p = rw.base_params(2)

    #define wrefernce frames
    if rot_AoA:
        fwt_root_frame = mf.HomogenousTransform().Translate(0,0,p.q[0]).R_y(p.alpha_r).R_x(-p.q[-1])
    else:
        fwt_root_frame = mf.HomogenousTransform().Translate(0,0,p.q[0]).R_x(-p.q[-1])
    fwt_flexural_frame = fwt_root_frame.Translate(p.x_f1,p.y_1,0)
    fwt_com_frame = fwt_root_frame.Translate(p.c/2,p.s_1/2,0)
    
    #Create Elemnts
    M_fwt = ele.MassMatrix(p.m_1,I_xx = p.I_xx_1)
    fwt_ele = ele.RigidElement(fwt_com_frame,M_fwt,True)

    spring_ele = ele.Spring(p.q[0],p.EI)
  
    alpha_fwt =  0
    alphadot_fwt = 0

    if fwt_free:
        alpha_fwt +=  p.alpha_1
        alphadot_fwt += p.alphadot_1

    if fwt_frot:
        fwt_AeroForces_perUnit = ef.AeroForce_2.PerUnitSpan(p,fwt_flexural_frame,p.a_1,
                            alphadot = p.alphadot_1,
                            M_thetadot = p.M_thetadot,
                            e = p.e_1,
                            rootAlpha = p.alpha_1,
                            alpha_zero = 0)
    else:
        fwt_AeroForces_perUnit = ef.AeroForce_1.PerUnitSpan(p,fwt_flexural_frame,p.a_1,
                            alphadot = p.alphadot_1,
                            M_thetadot = p.M_thetadot,
                            e = p.e_1,
                            rootAlpha = p.alpha_1,
                            deltaAlpha = 0, 
                            alpha_zero = 0)

    forces = []
    segments = 5
    for i in range(segments):
        seg_width = p.s_1/segments
        yi = seg_width/2 + i*seg_width
        forces.append(fwt_AeroForces_perUnit.subs({p.y_1:yi})*seg_width)
    Q = sym.Matrix([0]*p.qs)
    for f in forces:
        Q += f.Q()
    fwt_AeroForces = ef.ExternalForce(Q)

    # Setup AoA of FWT
    fwt_aoa = mf.GetAoA(p.alpha_r,0,p.Lambda,0 if not fwt_free else p.q[-1])

    ## Sub in Aero Forces
    fwt_AeroForces = fwt_AeroForces.subs({p.alpha_1:fwt_aoa,p.alphadot_1:fwt_aoa.diff(t)})

    # Create the SYmbolic Model
    sm = mf.SymbolicModel.FromElementsAndForces(p,[fwt_ele,spring_ele],fwt_AeroForces)

    return sm,p


