"""Unit tests for the PetersFinite aero model (finite-state inflow, Stage 2).

Each test is framed by the bug class it catches. Two layers:

  (A) inflow core -- inflow_matrices coefficients and the Theodorsen C(k) it
      reproduces. Independent of the structural coupling.
  (B) descriptor blocks -- M_a / C_a / K_aero_lambda / aero_forcing(_vel) and
      the (4+N) eigenproblem they assemble, against the QS circulatory baseline.

Sign conventions under test (see peters_finite.py docstring for the derivation):
    w   = U* alpha + xi' + (1/2-a) alpha'              (3/4-chord downwash)
    aero_forcing      S        = [+1, (1/2-a)]         (q'' -> inflow, E side)
    aero_forcing_vel  [0,-U*]                          (q'  -> inflow, A side)
    +lambda_0 in alpha_eff, lambda_0 = 0.5 b_bar . Lambda
"""

import numpy as np
import pytest
from scipy.special import hankel2

from mflco.model.params import TypicalSectionParameters
from mflco.model.michigan_params import (
    calibrate_michigan, section_from_params, structural_zeta,
)
from mflco.aero.quasi_steady import QuasiSteady
from mflco.aero.peters_finite import PetersFinite, inflow_matrices
from mflco.model.analysis import (
    modal_analysis, descriptor_matrices, linearized_eigenvalues,
    undamped_natural_frequencies,
)

p = TypicalSectionParameters()        # Isogai Case A (default)
a = p.a
mu = p.mu


# =====================================================================
# (A) inflow core
# =====================================================================
@pytest.mark.parametrize("N, expected", [
    (2, [2.0, -1.0]),
    (3, [6.0, -6.0, 1.0]),
    (5, [20.0, -90.0, 140.0, -70.0, 1.0]),
])
def test_b_bar_known_values(N, expected):
    """b_bar matches Peters' closed form. Catches an off-by-one in the
    factorial recurrence (the load-bearing weights for lambda_0)."""
    _, b_bar, _, _ = inflow_matrices(N)
    assert np.allclose(b_bar, expected)


def test_c_bar_is_two_over_n():
    """c_bar[n] = 2/n. Catches a wrong forcing-weight index."""
    _, _, c_bar, _ = inflow_matrices(4)
    assert np.allclose(c_bar, [2.0, 1.0, 2.0 / 3.0, 0.5])


def test_lag_poles_stable():
    """eig(A_bar^-1) all have Re>0, so the lag eigenvalues -u* eig(A_bar^-1) are
    stable (decaying). A_bar is not symmetric, so some poles are complex -- that
    is expected, not a bug; this only asserts decay."""
    for N in (3, 5, 8):
        A_bar, *_ = inflow_matrices(N)
        assert np.all(np.linalg.eigvals(np.linalg.inv(A_bar)).real > 0)


def test_high_freq_limit_half():
    """b_bar^T A_bar^-1 c_bar -> 1 as N grows, i.e. C_N(inf) -> 1/2 (Theodorsen's
    high-frequency limit). Catches a scaling error in A_bar/b_bar/c_bar."""
    A_bar, b_bar, c_bar, _ = inflow_matrices(12)
    val = b_bar @ np.linalg.solve(A_bar, c_bar)
    assert val == pytest.approx(1.0, abs=1e-3)


def _peters_Ck(N, k):
    """Finite-state lift-deficiency C_N(k) = 1 - 0.5 i k b^T (i k A_bar + I)^-1 c.
    Derived from the -wdot forcing with +lambda_0 feedback (the sign pairing that
    yields C<=1; the opposite pairing gives C>1, which is unphysical)."""
    A_bar, b_bar, c_bar, _ = inflow_matrices(N)
    M = 1j * k * A_bar + np.eye(N)
    return 1.0 - 0.5j * k * (b_bar @ np.linalg.solve(M, c_bar))


def _theodorsen(k):
    H1 = hankel2(1, k)
    H0 = hankel2(0, k)
    return H1 / (H1 + 1j * H0)


@pytest.mark.parametrize("k", [0.05, 0.1, 0.2, 0.4, 0.7, 1.0, 2.0, 5.0])
def test_Ck_recovers_theodorsen(k):
    """The whole point of Stage 2: N=8 inflow states reproduce Theodorsen's C(k)
    across reduced frequencies. Validates {c_bar, A_bar, b_bar} magnitude AND the
    forcing/feedback sign pairing together. Catches the qualitative failure that
    a flipped global sign would produce (C drifting to >1)."""
    cn = _peters_Ck(8, k)
    ct = _theodorsen(k)
    assert cn.real == pytest.approx(ct.real, abs=0.02)
    assert cn.imag == pytest.approx(ct.imag, abs=0.02)


def test_Ck_stays_physical():
    """|C_N(k)| <= ~1 everywhere (Theodorsen never amplifies). A flipped forcing
    sign sends C_N(k) above 1 -- this is the cheap guard against that."""
    for k in np.linspace(0.01, 10, 40):
        assert abs(_peters_Ck(8, k)) <= 1.01


# =====================================================================
# (B) descriptor blocks
# =====================================================================
def test_apparent_mass_symmetric_and_adds_inertia():
    """M_a symmetric; M_s - M_a is heavier on the diagonal than M_s (added mass
    increases both modal inertias). Catches a sign slip in M_a."""
    aero = PetersFinite(p, M_inf=0.0, N=3)
    Ma = aero.M_a()
    assert np.allclose(Ma, Ma.T)
    Meff = p.mass_matrix() - Ma
    assert Meff[0, 0] > p.mass_matrix()[0, 0]
    assert Meff[1, 1] > p.mass_matrix()[1, 1]


def test_apparent_mass_matches_standard_added_mass():
    """M_a = -(standard Theodorsen added mass [[1,-a],[-a,a^2+1/8]]/mu)."""
    aero = PetersFinite(p, M_inf=0.0, N=3)
    added = np.array([[1.0, -a], [-a, a ** 2 + 0.125]]) / mu
    assert np.allclose(aero.M_a(), -added)


def test_apparent_damping_sign():
    """REGRESSION: C_a (non-circulatory damping) must be POSITIVE on the alpha'
    column: C_a[0,1] = +U*/mu, C_a[1,1] = +U*(1/2-a)/mu. These come from the
    +U* alpha' term of Theodorsen's L_nc and the -U*(1/2-a) alpha' of M_nc, moved
    to the damping side. A negated C_a (the original draft) leaves the inflow C(k)
    untouched but quietly destabilises flutter (Michigan 13.1 -> 10.5 m/s), which
    masquerades as a missing-3D/viscous deficit. Pin the sign."""
    aero = PetersFinite(p, M_inf=0.0, N=3)
    u = 1.7
    Ca = aero.C_a(u)
    assert Ca[0, 1] == pytest.approx(u / mu)
    assert Ca[1, 1] == pytest.approx(u * (0.5 - a) / mu)
    assert Ca[0, 1] > 0 and Ca[1, 1] > 0          # the sign that matters
    assert np.allclose(Ca[:, 0], 0.0)             # no plunge-rate (xi') term


def test_apparent_terms_off_toggle():
    """include_apparent_mass=False zeroes both M_a and C_a (the U*=0 oracle path).
    Catches a toggle that forgets one of the two blocks."""
    aero = PetersFinite(p, M_inf=0.0, N=3, include_apparent_mass=False)
    assert np.allclose(aero.M_a(), 0.0)
    assert np.allclose(aero.C_a(2.0), 0.0)


def test_aero_forcing_sign_structure():
    """aero_forcing S = c_bar (x) [+1, (1/2-a)]: the plunge (xi'') and pitch
    (alpha'') entries share a sign (both trace to the same downwash w). The prior
    draft had them opposite. Also: the alpha'' entry of S and the alpha' entry of
    the velocity block are OPPOSITE (the cross-block tell of q'' moved across the
    descriptor '=')."""
    aero = PetersFinite(p, M_inf=0.0, N=3)
    c_bar = aero.c_bar
    F = aero.aero_forcing()                        # (N,2) = c_bar (x) S
    # recover S from the first inflow row (c_bar[0] = 2 != 0)
    S = F[0] / c_bar[0]
    assert S[0] == pytest.approx(1.0)
    assert S[1] == pytest.approx(0.5 - a)
    assert np.sign(S[0]) == np.sign(S[1])          # same sign (since 1/2-a>0 here)
    Fvel = aero.aero_forcing_vel(1.3)              # c_bar (x) [0,-U*]
    vel = Fvel[0] / c_bar[0]
    assert vel[0] == pytest.approx(0.0)
    assert vel[1] == pytest.approx(-1.3)
    assert np.sign(S[1]) == -np.sign(vel[1])       # opposite across the '='


def test_velocity_forcing_zero_at_zero_speed():
    """aero_forcing_vel vanishes at U*=0 -- which is exactly why the U*=0 oracle
    cannot test it (the Michigan sweep does). Documents the blind spot as a test."""
    aero = PetersFinite(p, M_inf=0.0, N=3)
    assert np.allclose(aero.aero_forcing_vel(0.0), 0.0)
    assert not np.allclose(aero.aero_forcing_vel(1.0), 0.0)


def test_circulatory_blocks_match_quasi_steady():
    """Peters' circulatory K_aero/C_aero are byte-for-byte the QS ones (Stage 2
    adds wake lag + apparent mass ON TOP of the same circulation, it does not
    re-derive it). Catches accidental divergence between the two models."""
    qs = QuasiSteady(p, M_inf=0.3)
    pet = PetersFinite(p, M_inf=0.3, N=4)
    for u in (0.5, 1.7, 3.0):
        assert np.allclose(qs.K_aero(u), pet.K_aero(u))
        assert np.allclose(qs.C_aero(u), pet.C_aero(u))


def test_K_aero_lambda_matches_forces():
    """K_aero_lambda is the inflow->structure force coupling in the descriptor A
    block; it must equal the inflow term of the time-marching force forces().
    With y_struct=0 the circulatory bracket is purely u* lambda_0, so column j of
    K_aero_lambda equals forces() evaluated at y_aero = e_j. Cross-checks Job A
    (descriptor coupling) against Job B (force model) and pins the sign that the
    C(k) feedback and Michigan flutter both depend on. (NB it is NOT the K_aero
    alpha-column scaled by 0.5 b_bar/u*: K_aero = -dQ/dalpha while K_aero_lambda =
    +dQ/dlambda, so that route carries an extra minus sign.)"""
    aero = PetersFinite(p, M_inf=0.0, N=3)
    u = 1.7
    Kl = aero.K_aero_lambda(u)                     # (2,N)
    zero_struct = [0.0, 0.0, 0.0, 0.0]
    for j in range(aero.N):
        ej = np.eye(aero.N)[j]
        Qj = aero.forces(0.0, zero_struct, ej, u)  # force from a unit inflow state
        assert np.allclose(Kl[:, j], Qj)


def test_descriptor_shape_and_nonsingular_E():
    """E, A are (4+N)^2 and E is invertible (det(M_s-M_a) det(A_bar) != 0), so the
    generalized eigenproblem A x = lambda E x is well posed."""
    aero = PetersFinite(p, M_inf=0.0, N=3)
    E, A = descriptor_matrices(p, 0.0, 1.5, aero)
    assert E.shape == (7, 7) and A.shape == (7, 7)
    assert abs(np.linalg.det(E)) > 1e-6


def test_eigenvalue_count():
    """(4+N) finite eigenvalues; exactly 2 structural oscillatory pairs are
    selectable, the rest are lag modes."""
    aero = PetersFinite(p, M_inf=0.0, N=3)
    ev = linearized_eigenvalues(p, 0.0, 1.5, aero)
    assert len(ev) == 7
    # 2 structural modes recoverable via modal_analysis (least-damped oscillatory)
    f, z = modal_analysis(p, 0.0, 1.5, aero)
    assert len(f) == 2 and len(z) == 2


# =====================================================================
# (B) U*=0 oracle on the calibrated Michigan section (dimensional check)
# =====================================================================
def test_michigan_u0_oracle_apparent_mass_off():
    """At U*=0 with apparent mass OFF, the inflow decouples (K_aero_lambda(0)=0)
    and the wind-off modes must equal the structural 5.3/6.2 Hz calibration."""
    cal = calibrate_michigan(zeta=structural_zeta())
    sec = section_from_params(cal)
    aero = PetersFinite(sec, M_inf=0.0, N=3, include_apparent_mass=False)
    f, _ = modal_analysis(sec, 0.0, 0.0, aero)
    hz = np.sort(f) * cal.omega_alpha / (2 * np.pi)
    assert hz[0] == pytest.approx(5.3, abs=0.02)
    assert hz[1] == pytest.approx(6.2, abs=0.02)


def test_apparent_mass_lowers_wind_off_frequencies():
    """Apparent mass adds inertia, so turning it ON must LOWER both wind-off
    frequencies relative to OFF. Catches an added-mass sign that raises them."""
    cal = calibrate_michigan(zeta=structural_zeta())
    sec = section_from_params(cal)
    f_off, _ = modal_analysis(sec, 0.0, 0.0,
                              PetersFinite(sec, 0.0, N=3, include_apparent_mass=False))
    f_on, _ = modal_analysis(sec, 0.0, 0.0,
                             PetersFinite(sec, 0.0, N=3, include_apparent_mass=True))
    assert np.sort(f_on)[0] <= np.sort(f_off)[0]
    assert np.sort(f_on)[1] <= np.sort(f_off)[1]