"""Peters' finite-state inflow model.

Essentially Quasi-Steady but with added wake-lag through the induced
velocity variable (lambda_zero), carried by N inflow states.

Convention matches quasi_steady.py: tau = omega_alpha*t, xi = h/b,
U* = U/(b*omega_alpha). Effective AoA gains +lambda_0/U* (FlowLab sign),
so the circulatory bracket gains one term:  circ += U* * lambda_0,
with lambda_0 = 0.5 * b_bar . Lambda.
"""

import numpy as np
from math import factorial


def inflow_matrices(N):
    """Peters finite-state inflow coefficient matrices (1-indexed n in theory).

    Returns
    -------
    A_bar : (N, N)  inflow ODE matrix (invertible; A_bar^-1 has Re>0 -> stable lag)
    b_bar : (N,)    induced-flow weights, lambda_0 = 0.5 * b_bar . Lambda
    c_bar : (N,)    forcing weights, c_bar[n] = 2/n
    d_bar : (N,)    [1/2, 0, 0, ...]
    Known: b_bar = [2,-1] (N=2), [6,-6,1] (N=3), [20,-90,140,-70,1] (N=5).
    """
    if N < 1:
        raise ValueError("N must be >= 1")
    n = np.arange(1, N + 1)
    c_bar = 2.0 / n
    d_bar = np.zeros(N)
    d_bar[0] = 0.5
    b_bar = np.empty(N)
    for k in range(1, N + 1):
        if k != N:
            b_bar[k - 1] = (-1) ** (k - 1) * factorial(N + k - 1) \
                / factorial(N - k - 1) / factorial(k) ** 2
        else:
            b_bar[k - 1] = (-1) ** (k - 1)
    D = np.zeros((N, N))
    for i in range(1, N + 1):
        for j in range(1, N + 1):
            if i == j + 1:
                D[i - 1, j - 1] = 1.0 / (2 * i)
            elif i == j - 1:
                D[i - 1, j - 1] = -1.0 / (2 * i)
    A_bar = (D + np.outer(d_bar, b_bar) + np.outer(c_bar, d_bar)
             + 0.5 * np.outer(c_bar, b_bar))
    return A_bar, b_bar, c_bar, d_bar


class PetersFinite:

    def __init__(self, params, M_inf, N=3, include_apparent_mass=True):
        self.params = params
        self.M_inf  = M_inf
        self.N      = N  # additional states; Peters' recommended range [3, 10]
        # toggle for the U*=0 oracle: M_a is U*-independent and shifts the modes
        # even at zero airspeed, so set False to recover the pure structural freqs
        self.include_apparent_mass = include_apparent_mass
        self.A_bar, self.b_bar, self.c_bar, self.d_bar = inflow_matrices(N)

    @property
    def n_aero_states(self) -> int:
        return self.N

    def K_aero(self, u_star):
        """(2,2) circulatory stiffness, identical to QuasiSteady.K_aero (Stage 2
        adds wake lag + apparent mass ON TOP of the same circulation)."""
        comp_factor = 1/np.sqrt(1 - self.M_inf**2)
        g = (2/self.params.mu)*comp_factor
        K_mat = np.asarray([[0, g*u_star**2],
                            [0, -g*(0.5 + self.params.a)*u_star**2]])
        return K_mat

    def C_aero(self, u_star):
        """(2,2) circulatory damping, identical to QuasiSteady.C_aero."""
        comp_factor = 1/np.sqrt(1 - self.M_inf**2)
        g = (2/self.params.mu)*comp_factor
        C_mat = np.asarray([[g*u_star,                       g*u_star*(0.5 - self.params.a)],
                            [-g*(0.5 + self.params.a)*u_star, -g*(0.5 + self.params.a)*u_star*(0.5 - self.params.a)]])
        return C_mat

    def K_aero_lambda(self, u_star):
        """(2, N) block: how the inflow states feed the structural forces.
        Matrix form of the 'circ += u_star*lambda_0' inflow term, with Lambda
        factored out.  K_aero_lambda @ Lambda gives the inflow generalized force."""
        mu, a = self.params.mu, self.params.a
        comp  = 1.0 / np.sqrt(1.0 - self.M_inf**2)
        row_xi    = (-2/mu)         * comp * u_star * 0.5 * self.b_bar   # (N,)
        row_alpha = (+2/mu)*(0.5+a) * comp * u_star * 0.5 * self.b_bar   # (N,)
        return np.array([row_xi, row_alpha])                            # (2, N)

    def M_a(self):
        """(2,2) apparent-mass matrix: coefficients of [xi'', alpha''] from the
        non-circulatory terms (sec 3.3). Symmetric (genuine fluid inertia).
        Constant: no U*, no PG (apparent mass is an inertial effect)."""
        if not self.include_apparent_mass:
            return np.zeros((2, 2))
        mu, a = self.params.mu, self.params.a
        return np.array([[-1/mu,        a/mu            ],
                         [ a/mu,   -(a**2 + 0.125)/mu   ]])

    def C_a(self, u_star):
        """(2,2) apparent-damping matrix: coefficients of [xi', alpha'] from the
        non-circulatory terms (sec 3.3). Only the alpha' column is nonzero
        (no xi' term in apparent mass). Scales with U*.

        SIGN: the entries are the alpha' coefficients moved onto the DAMPING side.
        forces()'s apparent_mass comment has +u_star*alpha' (xi-row) and
        -u_star*(0.5-a)*alpha' (alpha-row) on the FORCE side; moving those across
        to the C matrix flips them, giving +u_star/mu and +u_star*(0.5-a)/mu.
        (The force-side signs in the matrix destabilise Michigan ~13.1 -> 10.5.)"""
        if not self.include_apparent_mass:
            return np.zeros((2, 2))
        mu, a = self.params.mu, self.params.a
        return np.array([[0.0,   u_star/mu              ],
                         [0.0,   u_star*(0.5 - a)/mu    ]])

    def forces(self, tau, y_struct, y_aero, u_star):
        [xi, alpha, xi_dot, alpha_dot] = y_struct   # structural state

        mu = self.params.mu
        a  = self.params.a
        M_inf = self.M_inf
        b_bar = self.b_bar

        comp_factor = 1.0 / np.sqrt(1.0 - M_inf**2)
        circ_shared = (u_star**2 * alpha) + (u_star * xi_dot) \
            + (u_star * (0.5 - a) * alpha_dot)

        lambda_zero = 0.5 * (b_bar @ y_aero)        # y_aero IS Lambda; dot product
        induced_velocity = u_star * lambda_zero

        circ_total = circ_shared + induced_velocity  # inflow adds INSIDE the bracket
        Q_xi    = (-2/mu)         * circ_total * comp_factor
        Q_alpha = (+2/mu)*(0.5+a) * circ_total * comp_factor

        return np.array([Q_xi, Q_alpha])             # (2,)

        # apparent mass (Step 2 -> goes into M_a/C_a matrices, not here):
        # apparent_mass_xi    = xi_ddot + u_star*alpha_dot - a*alpha_ddot
        # apparent_mass_alpha = a*(xi_ddot - a*alpha_ddot) - u_star*(0.5-a)*alpha_dot - 0.125*alpha_ddot

    def aero_rhs(self, tau, y_struct, y_aero, u_star):
        """Inflow state derivatives Lambda' (N,) for the LCO time-march (Job B).

        Peters inflow ODE (FlowLab -wdot forcing):
            A_bar Lambda' = -u_star*Lambda - c_bar * w'
        where the 3/4-chord downwash rate w' = xi'' + u_star*alpha' + (0.5-a)*alpha''
        carries the structural ACCELERATIONS. Rather than thread q'' through this
        signature, recompute it here from the SAME effective EOM the solver uses
        (so the two stay consistent):
            (M_s - M_a) q'' + (C_s + C_a) q' + K_s q = forces()
        For the U*=0 / apparent-mass-off case M_a=C_a=0, so this is the bare
        structural solve. The 2x2 solve is cheap; keeping it self-contained also
        lets the solver's setup-time probe call aero_rhs before any q'' exists."""
        [xi, alpha, xi_dot, alpha_dot] = y_struct
        a = self.params.a
        q     = np.array([xi, alpha])
        q_dot = np.array([xi_dot, alpha_dot])

        M_eff = self.params.mass_matrix()    - self.M_a()         # (M_s - M_a)
        C_eff = self.params.damping_matrix() + self.C_a(u_star)   # (C_s + C_a)
        K     = self.params.stiffness_matrix(alpha)               # nonlinear in alpha
        Q     = self.forces(tau, y_struct, y_aero, u_star)        # circulatory (incl. lambda_0)
        q_ddot = np.linalg.solve(M_eff, Q - (C_eff @ q_dot) - (K @ q))

        # w' = S . q'' + u_star*alpha' , S = [1, 0.5-a]
        w_dot = q_ddot[0] + (0.5 - a)*q_ddot[1] + u_star*alpha_dot
        Lambda = np.asarray(y_aero, dtype=float)
        Lambda_dot = np.linalg.solve(self.A_bar, -u_star*Lambda - self.c_bar*w_dot)
        return Lambda_dot

    def aero_forcing(self):
        """(N, 2) block: how structural accelerations [xi'', alpha''] force the
        lambda-ODE.  c_bar (outer) S, with S = [+1, (0.5 - a)].  Goes on the E side.
        Both entries share a sign (both are components of the same 3/4-chord
        downwash w' = xi'' + u_star*alpha' + (0.5-a)*alpha''); the draft's [-1, ..]
        flipped the xi'' entry against the alpha'' entry and against forces()."""
        a = self.params.a
        S = np.array([1.0, 0.5 - a])      # acceleration-forcing row, (2,)
        return np.outer(self.c_bar, S)      # (N, 2)

    def aero_forcing_vel(self, u_star):
        """(N, 2) block: how structural velocities [xi', alpha'] force the
        lambda-ODE.  c_bar (outer) [0, -u_star] -- the u_star*alpha' term of w'
        (the rate of the steady-AoA part of the downwash).  Goes on the A side.
        Zero at U*=0, so the U*=0 oracle cannot test it (the sweep does)."""
        vel = np.array([0.0, -u_star])      # velocity-forcing row, (2,)
        return np.outer(self.c_bar, vel)    # (N, 2)