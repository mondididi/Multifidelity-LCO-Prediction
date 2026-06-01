"""Time integration wrapper around scipy.integrate.solve_ivp for 2-DOF pitch-plunge typical section."""    

import numpy as np
from scipy.integrate import solve_ivp
from .eom import structural_rhs
from ..aero.none import NoAero

def integrate(params, y0, tau_span, aero=None, U_star = None,
              method='RK45', rtol=1e-8, atol=1e-10, t_eval=None):
    
    """Integrate the 2-DOF typical section in nondimensional time.

    Thin wrapper around scipy.integrate.solve_ivp.

    Parameters
    ----------
    params : TypicalSectionParameters
        Section parameters; provides M, C, K matrices.
    y0 : array-like of length (4,) with aero ICs filled internally
        Initial state [xi, alpha, xi_dot, alpha_dot].
    tau_span : tuple of (float, float)
        Nondimensional time interval (tau_start, tau_end).
    aero : AeroModel object, not callable, (AeroModel or None)
        None defaults to NoAero (0 forces, 0 states) - stage 0 path
    U_star : Nondimensional airspeed - reduced velocity: U_star = U / (b*omega_alpha)
        bifurcation variable, sweep upward to watch for flutter boundary        
    method : str, default "RK45"
        solve_ivp integration method.
    rtol, atol : float
        Relative and absolute tolerances.
    t_eval : array-like or None, default None
        Times to store the solution. None lets the solver choose.

    Returns
    -------
    sol : OdeResult
        sol.t : (N,) nondimensional times.
        sol.y : (4+n_aero, N) state history — rows are xi, alpha, xi_dot, alpha_dot, [aero].
    """

    if aero is None:  #check for NoAero class
        aero = NoAero()

    #build augmented initial state
    y0_struct = np.asarray(y0, float) #(4,)
    y0_aero = np.zeros(aero.n_aero_states) #(n,): represents the wake's memory of the airfoil's past motion (lag variable)

    #single setup-time assertion
    probe = np.asarray(aero.aero_rhs(tau_span[0], y0_struct, y0_aero, U_star)) #initial state of rhs
    assert probe.shape == (aero.n_aero_states,), f"aero_rhs len {probe.shape} != n_aero_states {aero.n_aero_states}"
    #check is the shape of n_aero_states and aero_rhs agrees

    y0_full = np.concatenate([y0_struct,y0_aero]) #combine into augmented initial state

    sol = solve_ivp(
        structural_rhs,           
        tau_span,
        y0_full,
        method=method,
        rtol=rtol,
        atol=atol,
        t_eval=t_eval,
        args= (params, aero, U_star), #account for aero stuffs
        )

    return sol #(4+n, N)
    