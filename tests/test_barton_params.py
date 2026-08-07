"""
test_barton_params.py
=====================
Gate the Barton (AeroelasticCBC.jl model.jl) mapping: stage-0 machinery
verification with ZERO free parameters. Every input is published; the pipeline
wind-off modes must match the raw dimensional eigenproblem to machine
precision, the closed-form mapping identities must hold exactly, and the
linear flutter must land on the verified stage-0 values.

The slow fold known answer (Peters fold 16.92 vs full-Sears port 16.81 vs rig
~16) lives in examples/stage0_barton_test.py -- down-sweeps are too slow for
unit tests; here only the fast layers are load-bearing in CI.
"""

import numpy as np
import pytest

from mflco.model.barton_params import BartonCal, BASE, CASES
from mflco.model.analysis import modal_analysis, undamped_natural_frequencies
from mflco.aero.quasi_steady import QuasiSteady
from mflco.aero.peters_finite import PetersFinite


def _flutter_ms(cal, aero, us_lo, us_hi, n):
    """Linear flutter [m/s] from the min-damping zero crossing, narrow bracket."""
    p = aero.params
    prev, Uf = None, None
    for u in np.linspace(us_lo, us_hi, n):
        _, d = modal_analysis(p, 0.0, u, aero)
        m = float(np.min(d))
        if prev is not None and prev > 0.0 >= m:
            Uf = u
            break
        prev = m
    assert Uf is not None, "no flutter crossing inside the bracket"
    return float(cal.ustar_to_ms(Uf))


@pytest.mark.parametrize("case", [1, 2])
def test_windoff_matches_dimensional(case):
    """Pipeline U=0 modes equal the RAW dimensional eigenproblem (computed with
    no shared mflco code). Machine precision -- the zero-assumption gate on the
    container, M/K assembly and scaling. Tolerance 1e-9 Hz: observed 0 and
    8.9e-16 on two machines; 1e-9 leaves headroom for BLAS ordering."""
    cal = BartonCal(case)
    f_dim = cal.windoff_dimensional_hz()
    nat, _ = undamped_natural_frequencies(cal.section())
    f_pipe = np.sort(np.asarray(nat)) * cal.omega_alpha / (2 * np.pi)
    assert f_pipe == pytest.approx(f_dim, abs=1e-9)


@pytest.mark.parametrize("case", [1, 2])
def test_mapping_identities(case):
    """Every nondimensional group equals its closed form from the published
    dimensional numbers -- the mapping is derivation, not fitting."""
    P = dict(BASE, **CASES[case])
    cal = BartonCal(case)
    assert cal.beta == pytest.approx(P["k_a3"] / P["k_a"], rel=1e-14)
    assert cal.delta == pytest.approx(P["k_a2"] / (3 * P["k_a3"]), rel=1e-14)
    assert cal.mu == pytest.approx(P["m_T"] / (np.pi * P["rho"] * P["b"] ** 2),
                                   rel=1e-14)
    assert cal.x_alpha == pytest.approx(P["m_w"] * P["x_a"] / P["m_T"],
                                        rel=1e-14)
    assert cal.r_alpha_sq == pytest.approx(P["I_a"] / (P["m_T"] * P["b"] ** 2),
                                           rel=1e-14)
    # omega_ratio^2 = (k_h/m_T)/(k_a/I_a): both spring rates, one identity
    assert cal.omega_ratio ** 2 == pytest.approx(
        P["k_h"] * P["I_a"] / (P["m_T"] * P["k_a"]), rel=1e-14)


@pytest.mark.parametrize("case", [1, 2])
def test_delta_form_spans_bartons_free_quadratic(case):
    """Barton's FREE k_a2 alpha^2 is exactly the constrained trim form of
    params.py: the nondimensional restoring moment K[1,1]*alpha must expand as
    r2*(alpha + (k_a2/k_a) alpha^2 + (k_a3/k_a) alpha^3). Checked pointwise --
    this is the executable version of the 'delta spans Barton' claim."""
    P = dict(BASE, **CASES[case])
    cal = BartonCal(case)
    sec = cal.section()
    for a in np.radians([-12.0, -4.0, 3.0, 9.0, 15.0]):
        moment = sec.stiffness_matrix(a)[1, 1] * a
        expect = cal.r_alpha_sq * (a + (P["k_a2"] / P["k_a"]) * a ** 2
                                   + (P["k_a3"] / P["k_a"]) * a ** 3)
        assert moment == pytest.approx(expect, rel=1e-12)


def test_case2_qs_flutter_regression():
    """QS + Barton case-2 structure flutters at 23.22 m/s (stage-0 verified,
    both machines). Bracket step ~0.03 m/s; tolerance 0.15 covers grid."""
    cal = BartonCal(2)
    aero = QuasiSteady(cal.section(), 0.0)
    Uf = _flutter_ms(cal, aero, us_lo=7.8, us_hi=8.8, n=100)
    assert Uf == pytest.approx(23.22, abs=0.15)


def test_case2_peters_flutter_regression():
    """Peters N=6 + Barton case-2 structure flutters at 23.98 m/s -- 0.08% off
    the independent full-Sears port's 24.0: the headline stage-0 anchor."""
    cal = BartonCal(2)
    aero = PetersFinite(cal.section(), 0.0, N=6)
    Uf = _flutter_ms(cal, aero, us_lo=8.2, us_hi=8.9, n=70)
    assert Uf == pytest.approx(23.98, abs=0.15)


#verified:
# - wind-off pipeline modes == raw dimensional eigenproblem to <1e-9 Hz, both cases
# - all nondimensional groups equal their closed forms from published numbers
# - params.py's constrained delta form reproduces Barton's free quadratic exactly
# - QS flutter 23.22 and Peters flutter 23.98 m/s locked as regression anchors
#   (fold known answer 16.92 vs 16.81 gated by examples/stage0_barton_test.py)