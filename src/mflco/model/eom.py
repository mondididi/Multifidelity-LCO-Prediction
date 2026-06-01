"""Equations of motion for 2-DOF pitch-plunge typical section.

provides RHS function for scipy.solve_ivp (expects first order)

state vector y = [xi, alpha, xi_dot, alpha_dot,[+n_aero_states]]
xi = h/b (nondimensional plunge)
alpha = pitch angle (rad)
dots = derivative w.r.t. non-dimensional time tau = t * omega_alpha
"""

import numpy as np

def structural_rhs(tau, y, params, aero, U_star):
    """RHS side dy/tau for 2-DOF typical

    Parameters:
        tau : float, non-dimensional time (tau = t * omega_alpha)
        y = [xi, alpha, xi_dot, alpha_dot,[aero states]] : state vector (4+n_aero,)
        params : TypicalSectionParameters, contains system parameters and functions for M, C, K
        aero : AeroModel object, not callable
        U_star : Nondimensional airspeed - reduced velocity: U_star = U / (b*omega_alpha)

    Returns:
        dydtau : np.array, shape (4+n_aero_states,), time derivative of state vector w.r.t. tau
    """
    #unpack state vector into structures and aero
    y_struct = y[:4]    #struct first 4 terms
    y_aero = y[4:]      #aero after 4 terms, how many 

    Q = np.asarray(aero.forces(tau,y_struct, y_aero, U_star)) #aero force, (2,) output where Q = Q_xi, Q_alpha, if not provided, will be NoAero class, so return zeroes

    q = np.array([y_struct[0],y_struct[1]]) # xi, alpha   
    q_dot = np.array([y_struct[2],y_struct[3]]) #xi_dot, alpha_dot 1st order derivative

    M = params.mass_matrix()    #mass function, 0 arg
    C = params.damping_matrix() #damping function, 0 arg
    K = params.stiffness_matrix(y_struct[1]) #stiffness function, 1 arg, pitch angle (alpha), struct 2nd term

    rhs = Q - (C@q_dot) - (K@q)  #rhs of eom, where Mq_ddot = Q - Cq_dot - Kq (@ for matrix multiplication)
    q_ddot = np.linalg.solve(M, rhs)  #q double dot, 2nd order, solved from M(q_ddot) = rhs ; q_ddot s.t. Mq_ddot = rhs => q_ddot = M^-1 rhs, solved via np.linalg.solve for numerical stability

    struct_deriv = np.asarray([q_dot[0], q_dot[1], q_ddot[0], q_ddot[1]]) #xi_dot,alpha_dot,xi_ddot,alpha_ddot
    aero_deriv = np.asarray(aero.aero_rhs(tau, y_struct, y_aero, U_star)) #(n_aero_states) 

    return np.concatenate([struct_deriv, aero_deriv]) #return state vector with aero added