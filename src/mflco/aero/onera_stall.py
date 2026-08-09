"""ONERA dynamic-stall model (McAlister, Lambert & Petot, NASA TP 2399, 1984).

The paper splits every aerodynamic load F into two components:

    F1' + lambda*F1 = lambda*F_l + (lambda*s + sigma)*alpha' + s*alpha''   (2)
    F2'' + a*F2' + r*F2 = -(r*Delta + e*Delta')                            (3)

with F = F1 + F2, primes taken w.r.t. REDUCED (semi-chord) time, F_l the linear
extrapolation of the static load curve, F_s the static load, and

    Delta = F_l - F_s        the STALL DEFICIT (zero in attached flow).

Eq. (2) is first-order with a real negative pole -- the attached-flow response.
Eq. (3) is second-order with complex-conjugate poles -- the stalled response,
and it is what supplies HYSTERESIS: lift on the upstroke differs from lift on
the downstroke at the same incidence.

Implementation choice: eq. (2) is NOT re-implemented here. In attached flow the
ONERA model is equivalent to Theodorsen unsteady aerodynamics (Tang & Dowell,
via Lee, Price & Wong 1999 sec. 3.2.2), and PetersFinite is already a validated
finite-state realisation of exactly that physics (Michigan flutter -0.32%, sec.
4.2). This rung therefore = Peters (attached) + eq. (3) (separation), so the
ONLY physics added relative to the Peters rung is dynamic separation, and the
attached limit (Delta == 0) reduces to Peters IDENTICALLY -- gated in
tests/test_onera_stall.py. It also sidesteps the eq. (2) coefficients
(lambda = 0.25, s = 0.12 in the paper) whose remaining entries sigma, a, r, e
are published only as cubic splines in figures.

Time convention: the paper's reduced time is tau_s = U*t/b (semi-chord), while
this package integrates in tau = omega_alpha*t with U* = U/(b*omega_alpha), so
d/dtau_s = (1/U*) d/dtau. Eq. (3) in package time becomes

    F2'' = -a*U* F2' - r*U*^2 F2 - r*U*^2 Delta - e*U* Delta'

carried by two appended states [F2, F2'].

Stall deficit from the bundled polar, consistent with QSStall's trim
convention (perturbation measured from the trim point):

    Delta(alpha_eff) = 2*pi*alpha_eff - [Cl(trim + alpha_eff) - Cl(trim)]

so in steady state F2 -> -Delta and the total lift returns the static polar.

COEFFICIENT LAWS (Petot form): the coefficients are NOT constants -- they vary
with the instantaneous ABSOLUTE stall deficit Delta_abs = 2*pi*alpha - Cl(alpha)
(Sanchez Martinez 2018 sec. 2.2.4.4: "they can vary with the angle of attack";
McAlister TP-2399 fig. 16-19 splines). The standard published parameterisation:

    a = A0 + A2*Delta_abs^2
    r = (R0 + R2*Delta_abs^2)^2
    e = E2*Delta_abs^2

so in ATTACHED flow (Delta_abs = 0) the stalled equation is unforced (e = 0,
deficit zero) and the model is exactly the attached one -- the property the
constant-e version lacked. At this rig's trim the airfoil sits AT Cl_max with
Delta_abs(10 deg) ~ 0.34, so the laws place the stalled mode in its
mildly-stalled regime even for small oscillations, which is the physically
correct operating point.

PROVENANCE, stated honestly: A0 = 0.25 equals McAlister's linear-range lambda
(extracted from TP-2399 text); R0 = 0.20, A2 = 0.10, R2 = 0.10   # NASA TM-88917 App. B (Petot 1984), E2 = -0.07  # NEGATIVE per -(r*D + e*Ddot); TM-88917 are
the commonly quoted NACA 0012 constants of the Petot-form laws, not re-derived
here for the NACA 0020 -- the standing limitation of the rung. Constant-
coefficient overrides are kept for sensitivity studies (a 2^3-corner sweep of
the constant model showed outcomes ranging over LCO / none / blowup, which is
itself the argument for the laws). See examples/stage4_onera.py.

F_l convention: the deficit driving the forcing is measured against the SAME
2*pi circulation the attached (Peters) part produces, so the steady-state total
lift returns the measured static polar exactly; ONERA's textbook F_l (zero-
alpha polar-slope extrapolation) differs by (2*pi - Cl'(0))*alpha, absorbed
into the deficit definition without loss.
"""

import numpy as np
from scipy.interpolate import CubicSpline

from mflco.aero.peters_finite import PetersFinite
from mflco.aero.qs_stall import load_polar

# --- Petot-form law constants (NACA 0012 heritage; see docstring) ----------------
A0, A2   = 0.25, 0.10   # a = A0 + A2*D2 : stalled-mode damping vs deficit^2
R0, R2   = 0.20, 0.10  # TM-88917 App. B
E2       = -0.07  # NEGATIVE per -(r*D + e*Ddot); TM-88917         # e = E2*D2 : deficit-rate gain, ZERO in attached flow
A_STALL  = None         # constant-coefficient overrides for sensitivity runs;
R_STALL  = None         # None (default) -> use the Delta-dependent laws
E_STALL  = None
TRIM_DEG = 10.0         # rig mean AoA -- where the polar is sampled (Garcia Perez)


class ONERAStall(PetersFinite):
    """Peters finite-state inflow (attached) + ONERA eq. (3) (dynamic stall).

    State layout: y_aero = [Lambda (N inflow states), F2, F2_dot].
    Setting a_stall/r_stall/e_stall does not change the attached response; a
    linear polar (Cl = 2*pi*alpha) makes Delta vanish identically and recovers
    PetersFinite exactly.
    """

    def __init__(self, params, M_inf, N=6, alpha_deg=None, cl=None,
                 trim_deg=TRIM_DEG, a_stall=A_STALL, r_stall=R_STALL,
                 e_stall=E_STALL, include_apparent_mass=True):
        super().__init__(params, M_inf, N=N,
                         include_apparent_mass=include_apparent_mass)
        if alpha_deg is None:
            alpha_deg, cl = load_polar()
        # cubic-spline polar: np.interp's piecewise-linear kinks put slope
        # discontinuities into Delta' and stall the stiff integrator (~30x);
        # the spline gives a C2 deficit and an exact derivative
        self._ag       = np.radians(alpha_deg)          # table grid [rad]
        self._cl_spl   = CubicSpline(self._ag, np.asarray(cl, float))
        self._trim     = np.radians(trim_deg)
        self._cl_trim  = float(self._cl_spl(self._trim))
        self.a_stall   = a_stall            # None -> Petot-form laws
        self.r_stall   = r_stall
        self.e_stall   = e_stall
        # absolute deficit at trim: 2*pi*a_t - Cl(a_t); ~0.34 here (trim AT Cl_max)
        self._delta_trim = float(2.0 * np.pi * self._trim
                                 - self._cl_spl(self._trim))

    @property
    def n_aero_states(self) -> int:
        return self.N + 2                                # inflow + [F2, F2_dot]

    # --- stall deficit ------------------------------------------------------------
    def _dcl(self, alpha_eff):
        """Static lift increment from trim (spline-evaluated)."""
        return float(self._cl_spl(self._trim + alpha_eff)) - self._cl_trim

    def _delta(self, alpha_eff):
        """Stall deficit Delta = F_l - F_s: linear extrapolation minus static.

        Positive when the static curve has fallen below the attached line, i.e.
        exactly the lift the airfoil LOSES to separation at that incidence.
        Identically zero for a linear (2*pi) polar -> attached limit."""
        return 2.0 * np.pi * alpha_eff - self._dcl(alpha_eff)

    def _ddelta_dalpha(self, alpha_eff):
        """d(Delta)/d(alpha_eff), exact from the spline derivative."""
        return 2.0 * np.pi - float(self._cl_spl(self._trim + alpha_eff, 1))

    def _alpha_eff(self, y_struct, u_star):
        """Quasi-steady effective incidence driving the deficit.

        Same bracket as QuasiSteady/QSStall (circ = U*^2 * alpha_eff). The
        inflow states are deliberately NOT included here: Delta is a property
        of the sectional incidence seen by the boundary layer, and the ONERA
        coefficients were identified from prescribed-incidence oscillations."""
        xi, alpha, xi_dot, alpha_dot = y_struct
        return alpha + xi_dot / u_star + (0.5 - self.params.a) * alpha_dot / u_star

    def _law_coeffs(self, alpha_eff):
        """(a, r, e) of eq. (3): published NACA 0012 laws (Petot 1984;
        NASA TM-88917 App. B) on the ABSOLUTE deficit D = delta + delta_trim,
        a = A0 + A2 D^2, sqrt(r) = R0 + R2 D^2, e = E2 D^2, unless a constant
        override was set (constants mode, kept for the sensitivity study)."""
        d2 = (self._delta(alpha_eff) + self._delta_trim) ** 2
        a_c = self.a_stall if self.a_stall is not None else A0 + A2 * d2
        r_c = self.r_stall if self.r_stall is not None else (R0 + R2 * d2) ** 2
        e_c = self.e_stall if self.e_stall is not None else E2 * d2
        return a_c, r_c, e_c

    # --- forces -------------------------------------------------------------------
    def forces(self, tau, y_struct, y_aero, u_star):
        """Peters circulatory force + the stalled-load contribution F2.

        F2 is a lift-coefficient increment, so it scales exactly as QSStall's
        table lift: lift = U*^2/(pi*mu) * Cl * comp.  Sign follows the same
        convention as the parent (Q_xi = -lift)."""
        Q = super().forces(tau, y_struct, y_aero[:self.N], u_star)
        F2   = y_aero[self.N] if len(y_aero) > self.N else 0.0
        comp = 1.0 / np.sqrt(1.0 - self.M_inf ** 2)
        lift = (u_star ** 2 / (np.pi * self.params.mu)) * F2 * comp
        return Q + np.asarray([-lift, (0.5 + self.params.a) * lift])

    # --- aero states ---------------------------------------------------------------
    def aero_rhs(self, tau, y_struct, y_aero, u_star):
        """[Lambda' (N), F2', F2''] -- inflow ODE plus eq. (3) appended.

        The parent's aero_rhs is NOT reused: it recomputes q'' from the
        effective EOM to build the downwash rate, and that solve must see the
        TOTAL force including the stalled load, or the inflow states would lag
        an attached-only structure. The inflow ODE below is therefore the
        parent's, with forces() supplying attached + stall.

        Delta' is taken as (dDelta/dalpha)*alpha_dot: the deficit rate is driven
        by the pitch rate, consistent with the pure-pitch oscillations from
        which the coefficients were identified, and it avoids threading
        structural accelerations through this signature."""
        xi, alpha, xi_dot, alpha_dot = y_struct
        a = self.params.a
        q_dot = np.array([xi_dot, alpha_dot])
        M_eff = self.params.mass_matrix()    - self.M_a()
        C_eff = self.params.damping_matrix() + self.C_a(u_star)
        K     = self.params.stiffness_matrix(alpha)
        Q     = self.forces(tau, y_struct, y_aero, u_star)     # attached + stall
        q_ddot = np.linalg.solve(M_eff, Q - (C_eff @ q_dot)
                                 - (K @ np.array([xi, alpha])))
        w_dot  = q_ddot[0] + (0.5 - a) * q_ddot[1] + u_star * alpha_dot
        Lambda = np.asarray(y_aero[:self.N], dtype=float)
        lam_dot = np.linalg.solve(self.A_bar,
                                  -u_star * Lambda - self.c_bar * w_dot)
        F2, F2_dot = y_aero[self.N], y_aero[self.N + 1]
        alpha_dot  = y_struct[3]
        a_eff      = self._alpha_eff(y_struct, u_star)
        delta      = self._delta(a_eff)                  # perturbation deficit
        delta_dot  = self._ddelta_dalpha(a_eff) * alpha_dot
        a_c, r_c, e_c = self._law_coeffs(a_eff)
        F2_ddot = (-a_c * u_star * F2_dot
                   - r_c * u_star ** 2 * F2
                   - r_c * u_star ** 2 * delta
                   - e_c * u_star * delta_dot)
        return np.concatenate([lam_dot, [F2_dot, F2_ddot]])
