"""QS + static stall table -- the first viscous rung, one term changed from QS.

Physics change, and ONLY this change, relative to QuasiSteady:
    linear circulatory lift   Cl = 2*pi*alpha_eff
    becomes                   Cl = Cl_table(alpha_trim + alpha_eff)
                                   - Cl_table(alpha_trim)
with alpha_eff = alpha + xi_dot/U* + (0.5-a)*alpha_dot/U* -- the identical
circulatory bracket (circ = U*^2 * alpha_eff). Everything else is inherited
unchanged: moment structure (zero circulatory moment at a=-0.5, lift through
the EA at quarter-chord), no wake states, compressibility factor. Any change
in behaviour is therefore attributable to static separation alone. The
stalled Cm break is NOT modelled (lift-only stall); that is ONERA's rung.

Bundled polar: NACA 0020 at Re = 1.6e5, built by thickness interpolation of
the Sheldahl & Klimas (SAND80-2114) NACA 0018 and 0021 tables as distributed
with Sandia's CACTUS (test/Airfoil_Section_Data/):
Cl_0020 = (1/3) Cl_0018 + (2/3) Cl_0021. Re matches the Michigan rig
(U ~ 12.5 m/s, c = 0.2 m -> Re ~ 1.7e5). Regeneration recipe lives in
examples/stage3_stall_probe.py.

Verification: with a pure-linear table (Cl = 2*pi*alpha) and zero trim,
forces() reproduces QuasiSteady to machine precision (observed 4e-17) --
gated in tests/test_qs_stall.py.
"""

import os

import numpy as np

from mflco.aero.quasi_steady import QuasiSteady

DEFAULT_POLAR = os.path.join(os.path.dirname(__file__), "data",
                             "naca0020_re1p6e5.csv")


def load_polar(path=DEFAULT_POLAR):
    """Return (alpha_deg, Cl) arrays from a provenance-headed polar CSV."""
    d = np.loadtxt(path, delimiter=",", comments="#")
    return d[:, 0], d[:, 1]


class QSStall(QuasiSteady):
    """QuasiSteady with the circulatory lift read from a static polar about
    a nonzero trim angle. trim_deg is where the polar is sampled -- for the
    Michigan rig this is the 10 deg mean AoA, which is where the trim angle
    finally enters the aerodynamics."""

    def __init__(self, params, M_inf, alpha_deg, cl, trim_deg):
        super().__init__(params, M_inf)
        self._ag = np.radians(alpha_deg)                  # table grid [rad]
        self._cl = cl
        self._trim = np.radians(trim_deg)
        self._cl_trim = float(np.interp(self._trim, self._ag, self._cl))
        # local slope at trim [1/rad], central difference on the table
        h = np.radians(0.25)
        self.slope_trim = float(
            (np.interp(self._trim + h, self._ag, self._cl)
             - np.interp(self._trim - h, self._ag, self._cl)) / (2 * h))

    def _dcl(self, alpha_eff):
        return float(np.interp(self._trim + alpha_eff, self._ag, self._cl)) \
            - self._cl_trim

    def forces(self, tau, y_struct, y_aero, U_star):
        xi, alpha, xi_dot, alpha_dot = y_struct
        mu, a = self.params.mu, self.params.a
        comp = 1.0 / np.sqrt(1.0 - self.M_inf ** 2)
        a_eff = alpha + xi_dot / U_star + (0.5 - a) * alpha_dot / U_star
        # linear check: _dcl -> 2*pi*a_eff makes lift = (2/mu)*U*^2*a_eff = QS
        lift = (U_star ** 2 / (np.pi * mu)) * self._dcl(a_eff) * comp
        return np.asarray([-lift, (0.5 + a) * lift])

    # linearisation about trim for modal_analysis: slope 2*pi -> slope_trim
    def K_aero(self, U_star):
        return (self.slope_trim / (2 * np.pi)) * super().K_aero(U_star)

    def C_aero(self, U_star):
        return (self.slope_trim / (2 * np.pi)) * super().C_aero(U_star)