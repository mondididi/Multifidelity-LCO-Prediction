"""
test_onera_stall.py
===================
Gate the ONERA rung construction: Peters (attached, validated) + the eq. (3)
stalled-load states must reduce to PetersFinite EXACTLY in the attached limit,
the stall deficit must vanish at trim and carry the polar's slope collapse,
the coefficient laws must reproduce the published NACA 0012 fit (Petot 1984;
NASA TM-88917 App. B), and the steady state of eq. (3) must cancel the
deficit so the total quasi-static lift returns the static polar.

The slow bifurcation results (bistable window, fold) live in
examples/stage4_onera.py -- kick probes and down-sweeps are too slow for unit
tests; here the fast construction layers are load-bearing in CI.
"""

import numpy as np
import pytest

from mflco.model.params import TypicalSectionParameters
from mflco.aero.peters_finite import PetersFinite
from mflco.aero.onera_stall import ONERAStall


def _section():
    """Any well-posed section: the construction identities are parameter-free."""
    return TypicalSectionParameters(a=-0.5, x_alpha=0.1, r_alpha_sq=0.5,
                                    omega_ratio=0.8, mu=150.0)


def _linear_pair(p):
    """(PetersFinite, ONERAStall with linear polar + zero trim) -- Delta == 0."""
    lin = np.linspace(-60, 60, 481)
    pet = PetersFinite(p, 0.0, N=6)
    on = ONERAStall(p, 0.0, N=6, alpha_deg=lin,
                    cl=2 * np.pi * np.radians(lin), trim_deg=0.0)
    return pet, on


def test_attached_limit_reduces_to_peters():
    """Linear table + zero trim makes the deficit vanish identically, so
    forces() and aero_rhs() must equal PetersFinite with F2 = F2' = 0.
    Tolerance 1e-13: observed 0.0 exactly -- the stall extension adds NOTHING
    in attached flow, which is what makes this a one-change rung."""
    p = _section()
    pet, on = _linear_pair(p)
    rng = np.random.default_rng(0)
    for _ in range(100):
        ys = rng.normal(0, 0.15, 4)
        lam = rng.normal(0, 0.05, 6)
        ya = np.concatenate([lam, [0.0, 0.0]])
        U = rng.uniform(1.0, 5.0)
        assert on.forces(0, ys, ya, U) == pytest.approx(
            pet.forces(0, ys, lam, U), abs=1e-13)
        assert on.aero_rhs(0, ys, ya, U)[:6] == pytest.approx(
            pet.aero_rhs(0, ys, lam, U), abs=1e-13)


def test_state_count_is_peters_plus_two():
    """y_aero = [Lambda (N), F2, F2'] -- the ONLY state-space change."""
    on = ONERAStall(_section(), 0.0, N=6)
    assert on.n_aero_states == 8


def test_deficit_zero_at_trim_with_collapsed_slope():
    """Delta(0) = 0 by construction (deficit measured FROM trim), and its
    slope there is 2*pi minus the polar's trim slope ~0.63 /rad -- the same
    90% collapse the stall rung found. Tolerance 0.05 on the slope: spline
    representation of the tabulated grid."""
    on = ONERAStall(_section(), 0.0, N=6)          # bundled 0020 polar, trim 10
    assert on._delta(0.0) == pytest.approx(0.0, abs=1e-12)
    assert on._ddelta_dalpha(0.0) == pytest.approx(2 * np.pi - 0.631, abs=0.05)
    assert on._delta(np.radians(5.0)) > 0.3        # deficit grows post-trim


def test_laws_match_published_naca0012_fit():
    """a = 0.25 + 0.10 D^2, sqrt(r) = 0.20 + 0.10 D^2, e = -0.07 D^2
    (Petot 1984; NASA TM-88917 App. B), evaluated on the ABSOLUTE deficit
    D = delta(alpha_eff) + delta_trim. Checked at alpha_eff = 0 against the
    closed form with D = delta_trim -- exact, the law is a formula."""
    on = ONERAStall(_section(), 0.0, N=6)
    D = on._delta(0.0) + on._delta_trim
    a, r, e = on._law_coeffs(0.0)
    assert a == pytest.approx(0.25 + 0.10 * D ** 2, rel=1e-12)
    assert r == pytest.approx((0.20 + 0.10 * D ** 2) ** 2, rel=1e-12)
    assert e == pytest.approx(-0.07 * D ** 2, rel=1e-12)


def test_steady_state_f2_cancels_deficit():
    """At frozen incidence (alpha_dot = 0) eq. (3) equilibrates at
    F2* = -Delta exactly: r*U*^2*F2 + r*U*^2*Delta = 0. So with
    F2 = -Delta and F2' = 0 the returned F2'' must vanish, meaning the
    quasi-static total lift falls back onto the static polar -- the model
    only ADDS physics dynamically."""
    on = ONERAStall(_section(), 0.0, N=6)
    for a_deg in (2.0, 6.0, 10.0):
        al = np.radians(a_deg)
        ys = np.array([0.0, al, 0.0, 0.0])
        ya = np.zeros(8)
        ya[6] = -on._delta(al)                     # F2 = -Delta(alpha_eff)
        F2_ddot = on.aero_rhs(0.0, ys, ya, u_star=3.0)[-1]
        assert F2_ddot == pytest.approx(0.0, abs=1e-12)


#verified:
# - attached limit (linear polar, zero trim) == PetersFinite exactly
#   (observed 0.0; the rung is a true one-change extension)
# - state layout N+2; deficit zero at trim with the 90% slope collapse
# - coefficient laws reproduce the published Petot/TM-88917 NACA 0012 fit
# - eq. (3) steady state cancels the deficit: statics return the polar
#   (bistable window / fold gated by examples/stage4_onera.py)
