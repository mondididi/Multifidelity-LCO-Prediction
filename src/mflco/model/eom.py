"""Equations of motion for 2-DOF pitch-plunge typical section.

provides RHS function for scipy.solve_ivp (expects first order)

state vector y = [xi, alpha, xi_dot, alpha_dot]
xi = h/b (nondimensional plunge)
alpha = pitch angle (rad)
dots = derivative w.r.t. non-dimensional time tau = t * omega_alpha
"""

import numpy as np

def structural_rhs(tau, y, params, aero_force=None):
    """RHS side dy/tau for 2-DOF typical

    Parameters:
        tau : float, non-dimensional time (tau = t * omega_alpha)
        y = [xi, alpha, xi_dot, alpha_dot] : state vector (4,)
        params : TypicalSectionParameters, contains system parameters and functions for M, C, K
        aero_force : callable, optional, function of (tau, y) that returns aerodynamic force

    Returns:
        dydtau : np.array, shape (4,), time derivative of state vector w.r.t. tau
    """

    xi, alpha, xi_dot, alpha_dot = y  # unpack state vector
    
    q = np.array([xi, alpha])   
    q_dot = np.array([xi_dot, alpha_dot]) #1st order derivative

    M = params.mass_matrix()    #mass function, 0 arg
    C = params.damping_matrix() #damping function, 0 arg
    K = params.stiffness_matrix(alpha) #stiffness function, 1 arg, pitch angle (alpha)

    if aero_force is None: #if aero_force is not provided, assume zero
        Q = np.zeros(len(q))  # zero aerodynamic force, same dimension as q [2x2 bc 2 DOF pitch plunge]
    else:  #calculate
        Q = np.asarray(aero_force(tau, y)) #convert force to np array for np manipulation || aero_force(tau,y) callable funct

    rhs = Q - (C@q_dot) - (K@q)  #rhs of eom, where Mq_ddot = Q - Cq_dot - Kq (@ for matrix multiplication)
    q_ddot = np.linalg.solve(M, rhs)  #q double dot, 2nd order, solved from M(q_ddot) = rhs ; q_ddot s.t. Mq_ddot = rhs => q_ddot = M^-1 rhs, solved via np.linalg.solve for numerical stability

    return np.array([xi_dot, alpha_dot, q_ddot[0], q_ddot[1]])  #return state derivative vector