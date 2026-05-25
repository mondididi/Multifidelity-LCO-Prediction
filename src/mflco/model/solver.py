"""Time integration wrapper around scipy.integrate.solve_ivp for 2-DOF pitch-plunge typical section."""    

import numpy as np
from scipy.integrate import solve_ivp
from .eom import structural_rhs

def integrate(params, y0, tau_span, aero_force=None,
              method='RK45', rtol=1e-8, atol=1e-10, t_eval=None):
    
    """Integrate the 2-DOF typical section in nondimensional time.

    Thin wrapper around scipy.integrate.solve_ivp.

    Parameters
    ----------
    params : TypicalSectionParameters
        Section parameters; provides M, C, K matrices.
    y0 : array-like of length 4
        Initial state [xi, alpha, xi_dot, alpha_dot].
    tau_span : tuple of (float, float)
        Nondimensional time interval (tau_start, tau_end).
    aero_force : callable or None, default None
        Function (tau, y) -> (Q_xi, Q_alpha). None for Stage 0 (no aero).
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
        sol.y : (4, N) state history — rows are xi, alpha, xi_dot, alpha_dot.
    """

    y0_array = np.asarray(y0, dtype=float) #convert to numpy float array

    sol = solve_ivp(
        structural_rhs,           
        tau_span,
        y0_array,
        method=method,
        rtol=rtol,
        atol=atol,
        t_eval=t_eval,
        args= (params, aero_force),
        )

    return sol
    