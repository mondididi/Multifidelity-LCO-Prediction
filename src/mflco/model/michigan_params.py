"""
michigan_params.py
==================
Structural calibration of the 2-DOF pitch-plunge section to the García Pérez
et al. (AIAA J) experimental rig -- "the Michigan paper".

Calibrates the *effective* structure to the measured U=0 modes (5.3 / 6.2 Hz):
geometry (a, x_alpha, mu) is held fixed; r_alpha_sq and omega_alpha are solved
so the section's U=0 coupled modes match the experiment. We never tune a
structural parameter to a U>0 quantity -- the aero-fidelity result lives in the
U>0 comparison.

Full reference (provenance, equations, derivation, sensitivity, test rationale):
    docs/michigan_params.md
"""

from __future__ import annotations
from dataclasses import dataclass
import numpy as np

# ----------------------------------------------------------------------------
# Verified rig numbers (García Pérez et al.; NACA 0020, Bristol Prandtl tunnel)
# ----------------------------------------------------------------------------
M_KG          = 2.4        # total wing mass                       [kg]
S_ALPHA_KGM   = 0.0285     # static moment about elastic axis      [kg.m]
IA_BARE_KGM2  = 3.401e-3   # bottom-up pitch inertia about EA      [kg.m^2]
KH_NOM_NM     = 1580.0     # nominal plunge stiffness (1.58 N/mm)  [N/m]
KA_NOM_NMRAD  = 11.53      # nominal pitch stiffness               [N.m/rad]
B_M           = 0.1        # semichord (chord 0.2 m)               [m]
SPAN_M        = 0.44       # span                                  [m]
RHO_AIR       = 1.225      # air density                           [kg/m^3]

# Measured U=0 coupled modes (Fig. 6 intercepts) -- the calibration targets
F_PLUNGE_HZ   = 5.3        # lower (plunge-led) branch at U=0       [Hz]
F_PITCH_HZ    = 6.2        # upper (pitch-led)  branch at U=0       [Hz]

# Geometry-fixed nondimensionals (NOT calibrated)
A_NONDIM      = -0.5                       # spar at 25% chord, EA aft of midchord
X_ALPHA       = S_ALPHA_KGM / (M_KG * B_M) # static imbalance = 0.11875
MU            = M_KG / (np.pi * RHO_AIR * B_M**2 * SPAN_M)  # mass ratio ~= 141.7


@dataclass(frozen=True)
class MichiganParams:
    """Calibrated effective structural parameters + dimensional scales."""
    # nondimensional (feed these to the section)
    a: float
    x_alpha: float
    r_alpha_sq: float
    mu: float
    omega_ratio: float
    zeta: float
    # dimensional scales (for plotting / overlaying experimental data)
    b: float            # semichord [m]
    omega_alpha: float  # effective uncoupled pitch frequency [rad/s]
    # effective dimensional structure (for the record / sanity checks)
    Ia_eff: float       # [kg.m^2]
    Kh_eff: float       # [N/m]
    Ka_eff: float       # [N.m/rad]

    # ---- dimensional conversion helpers -------------------------------------
    def omega_to_hz(self, omega_over_omega_alpha):
        """nondim frequency (omega/omega_alpha) -> Hz."""
        return np.asarray(omega_over_omega_alpha) * self.omega_alpha / (2 * np.pi)

    def ustar_to_ms(self, u_star):
        """nondim velocity U* = U/(b omega_alpha) -> m/s."""
        return np.asarray(u_star) * self.b * self.omega_alpha

    def ms_to_ustar(self, u_ms):
        """m/s -> U*."""
        return np.asarray(u_ms) / (self.b * self.omega_alpha)


def calibrate_michigan(f_lo_hz: float = F_PLUNGE_HZ,
                       f_hi_hz: float = F_PITCH_HZ,
                       omega_ratio: float = 1.0,
                       zeta: float = 0.0) -> MichiganParams:
    """
    Closed-form structural calibration to the two measured U=0 modes.

    Holds geometry fixed (a, x_alpha, mu); solves r_alpha_sq (from the mode
    ratio) and omega_alpha (from the absolute scale) so the U=0 coupled modes
    equal (f_lo, f_hi). The one-parameter family is pinned by `omega_ratio`
    = omega_h/omega_alpha (default 1.0 = exact veering; ~0.95 makes the lower
    mode plunge-led, matching the paper's branch labelling).

    Derivation (eigenproblem, the rho->T->r_alpha_sq closed form, scale, and
    back-out): see docs/michigan_params.md section 5.
    """
    if f_hi_hz <= f_lo_hz:
        raise ValueError("f_hi_hz must exceed f_lo_hz")

    x  = X_ALPHA
    wr = float(omega_ratio)
    rho_ratio = (f_lo_hz / f_hi_hz) ** 2          # = lam_lo / lam_hi  (known)
    T = rho_ratio + 2.0 + 1.0 / rho_ratio         # = (lam_lo+lam_hi)^2/(lam_lo lam_hi)

    denom = T * wr**2 - (wr**2 + 1.0) ** 2
    if denom <= 0:
        raise ValueError(
            f"No positive r_alpha_sq for omega_ratio={wr} and this mode ratio "
            f"(modes too far apart for this pin). Choose omega_ratio closer to 1."
        )
    r_alpha_sq = T * wr**2 * x**2 / denom

    # solve the nondim quadratic for the two eigenvalues lam = (omega/omega_alpha)^2
    A = r_alpha_sq - x**2
    B = -r_alpha_sq * (wr**2 + 1.0)
    C = r_alpha_sq * wr**2
    lam = np.sort(np.roots([A, B, C]).real)        # [lam_lo, lam_hi]

    f_alpha_hz  = f_lo_hz / np.sqrt(lam[0])         # == f_hi_hz / sqrt(lam[1])
    omega_alpha = 2 * np.pi * f_alpha_hz

    # back out effective dimensional structure (record only)
    Ia_eff = r_alpha_sq * M_KG * B_M**2
    Ka_eff = Ia_eff * omega_alpha**2
    Kh_eff = M_KG * (wr * omega_alpha) ** 2

    return MichiganParams(
        a=A_NONDIM, x_alpha=x, r_alpha_sq=float(r_alpha_sq), mu=float(MU),
        omega_ratio=wr, zeta=float(zeta),
        b=B_M, omega_alpha=float(omega_alpha),
        Ia_eff=float(Ia_eff), Kh_eff=float(Kh_eff), Ka_eff=float(Ka_eff),
    )


def michigan_section(omega_ratio: float = 1.0, zeta: float = 0.0):
    """
    Build the calibrated Michigan typical section.

    NOTE: wire this to the real mflco section constructor. Replace the import
    and call below with your actual class (kept explicit so the mapping from
    nondim params -> section is reviewable, not hidden).
    """
    p = calibrate_michigan(omega_ratio=omega_ratio, zeta=zeta)

    # --- INTEGRATION POINT (adjust to your package) --------------------------
    # from mflco.structural import TypicalSection
    # return TypicalSection(a=p.a, x_alpha=p.x_alpha, r_alpha_sq=p.r_alpha_sq,
    #                       mu=p.mu, omega_ratio=p.omega_ratio, zeta=p.zeta)
    # ------------------------------------------------------------------------
    return p  # returns MichiganParams until wired; carries the dim scales too


if __name__ == "__main__":
    p = calibrate_michigan()
    print("Calibrated Michigan section (omega_ratio = 1.0, veering pin)")
    print(f"  a           = {p.a}")
    print(f"  x_alpha     = {p.x_alpha:.5f}   (geometry)")
    print(f"  r_alpha_sq  = {p.r_alpha_sq:.5f}   (effective; bottom-up was "
          f"{IA_BARE_KGM2/(M_KG*B_M**2):.4f})")
    print(f"  mu          = {p.mu:.3f}   (geometry)")
    print(f"  omega_ratio = {p.omega_ratio:.5f}")
    print(f"  zeta        = {p.zeta}   (set via recovery-rate back-out)")
    print(f"  --- scales ---")
    print(f"  b           = {p.b} m")
    print(f"  omega_alpha = {p.omega_alpha:.4f} rad/s ({p.omega_alpha/2/np.pi:.4f} Hz)")
    print(f"  U[m/s]      = {p.b*p.omega_alpha:.4f} * U*")
    print(f"  I_eff/I_bare= {p.Ia_eff/IA_BARE_KGM2:.2f}x ,  "
          f"Kh_eff/Kh_nom={p.Kh_eff/KH_NOM_NM:.2f}x ,  "
          f"Ka_eff/Ka_nom={p.Ka_eff/KA_NOM_NMRAD:.2f}x")
    print(f"  Fig17 fold 11.85 m/s -> U*={p.ms_to_ustar(11.85):.3f} ;  "
          f"flutter 13.19 m/s -> U*={p.ms_to_ustar(13.19):.3f}")