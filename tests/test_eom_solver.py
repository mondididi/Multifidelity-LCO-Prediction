"""Tests for eom.py - runs free-vibration time integrations and asserts energy is conserved.

with no aero, no damping, the system is a conservative oscillator. Therefore, E_total must be constant.
If wrong, (ex: sign error, asym M leak, wrong damp, int. tol. too lose, etc) E_total will change.
"""

import pytest #note: def test_***function_name***() is a pytest convention
from mflco.model.params import TypicalSectionParameters
from mflco.model.solver import integrate
import numpy as np


def _total_energy(y, params): #use _ to mark as internal
    """Helper function to compute total mechanical energy of a single y state vector"""
    [xi, alpha, xi_dot, alpha_dot] = y
    q_dot = np.array([xi_dot,alpha_dot])
    M = params.mass_matrix()
    KE = 0.5 * q_dot @ M @ q_dot        #refer to structural model section 5, symmetric lagrangian form
    PE = (0.5 * params.omega_ratio**2 * xi**2) + (0.5 * params.r_alpha_sq * alpha**2) + (0.25 * params.r_alpha_sq * params.beta * alpha**4)
    E_total = KE + PE
    return  E_total


def test_energy_conservation_linear_case():
    """with beta = 0, and damping = 0, total energy should be conserved during free vibration."""
    p = TypicalSectionParameters(beta=0.0, zeta_h = 0.0, zeta_alpha = 0.0) #called p as an object, with no beta & damping
    y0 = [0.0, 0.05, 0.0, 0.0]              #initial condition: y = [xi, alpha, xi_dot, alpha_dot], (plunge, pitch, plunge_dot, pitch_dot)
    tau_span =  (0.0, 20.0)                 #nondimensional time span for integration, from 0 to 20
    t_eval = np.linspace(tau_span[0], tau_span[1], 501)    #501 grid points for time grid
    sol = integrate(p, y0, tau_span,rtol=1e-10, atol=1e-12,t_eval=t_eval) #integrate the system with very tight tolerances to minimize numerical energy drift, output t,y
    
    E_history = []
    for i in range(sol.y.shape[1]): #loop over sol.y (which is 4 x 501), check energy at each time step
        y = sol.y[:, i] #current state vector at time t
        E_history.append(_total_energy(y, p)) #compute total energy at current state and append to E_history
        #output E = (1x501)
   
    E_history = np.asarray(E_history)
    drift = np.abs(E_history-E_history[0]).max() /E_history[0] #compute max relative energy drift
    assert drift < 1e-6 #check if abs error > 1e-6


def test_energy_conservation_nonlinear_case():
    """similar to test 1 but with cubic spring, so beta > 0."""
    p = TypicalSectionParameters(beta=3.0, zeta_h = 0.0, zeta_alpha = 0.0) #called p as an object, with beta > 0, let = 3.0
    y0 = [0.0, 0.1, 0.0, 0.0]              #add pitch magnitude so cubic has effects
    tau_span =  (0.0, 20.0)                 #nondimensional time span for integration, from 0 to 20
    t_eval = np.linspace(*tau_span, 501)    #501 grid points for time grid, use *tau_span for unpacking
    sol = integrate(p, y0, tau_span,rtol=1e-10, atol=1e-12,t_eval=t_eval) #integrate the system with very tight tolerances to minimize numerical energy drift, output t,y
    
    E_history = []
    for i in range(sol.y.shape[1]): #loop over sol.y (which is 4 x 501), check energy at each time step
        y = sol.y[:, i] #current state vector at time t
        E_history.append(_total_energy(y, p)) #compute total energy at current state and append to E_history
        #output E = (1x501)

    E_history = np.asarray(E_history)
    drift = np.abs(E_history-E_history[0]).max() /E_history[0] #compute max relative energy drift
    assert drift < 1e-6 #check if abs error > 1e-6
    

def test_undamped_amplitude_no_decay():
    """suppose amplitude drift weirdly. test checks if damping = 0, oscillation maintains amp."""
    #since we assume linear case, frequency is constant in time. (linear oscillator) if non-linear, frequecny depends on amplitude
    p = TypicalSectionParameters(beta=0.0, zeta_h = 0.0, zeta_alpha = 0.0) #default params, linear case
    y0 = [0.0, 0.05, 0.0, 0.0]              
    tau_span =  (0.0, 20.0)
    t_eval = np.linspace(*tau_span, 2001)    #denser time grid to better capture amplitude
    sol = integrate(p, y0, tau_span,rtol=1e-10, atol=1e-12,t_eval=t_eval)

    alpha_history = sol.y[1, :] #extract pitch history (alpha) from sol.y, only first row
    midpoint_id = len(alpha_history) // 2 #find midpoint index
    amp_first = np.abs(alpha_history[:midpoint_id]).max() #all entries upto midpoint, exclusive
    amp_second = np.abs(alpha_history[midpoint_id:]).max() #all entries from midpoint to end, inclusive

    #if energy is conserved and nothing damps oscillation, both halves should have same max amplitude
    assert amp_second == pytest.approx(amp_first, rel=1e-3) #check if second half amplitude is approx equal to first half, with relative tolerance of 1e-6


# Verified by this test file:
# - Mass matrix M is symmetric Lagrangian form (asymmetric form would leak energy).
# - EOM right-hand side in eom.py correctly couples M, C, K with state vector.
# - Time integration via solver.py preserves the conservative structure (no hidden numerical dissipation at rtol=1e-10, atol=1e-12).
# - Linear case (β=0): total energy E = KE + PE drifts < 1 part in 10⁶ over 20 nondimensional time units.
# - Nonlinear case (β=3): total energy conserved with cubic spring active.
#   Confirms K[1,1] = r_α²(1 + β α²) in params.py matches PE = ¼ r_α² β α⁴
#   in this file — both implementations of the cubic term agree.
# - With ζ_h = ζ_α = 0, oscillation amplitude does not decay over 20 time units — no hidden damping, no spurious energy growth.