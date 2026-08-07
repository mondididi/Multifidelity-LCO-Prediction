"""
test_qs_stall.py
================
Gate the stall rung construction: the one-term extension of QuasiSteady must
reduce to QS exactly under a linear table, and the bundled NACA 0020 polar
must carry the physical finding the rung's result rests on -- the rig trims
at its airfoil's static Cl_max, collapsing the local lift slope by ~90%.

The slow probe result (no linear flutter, no finite-amplitude LCO anywhere)
lives in examples/stage3_stall_probe.py -- kick grids are too slow for unit
tests; here the fast construction layers are load-bearing in CI.
"""

import numpy as np
import pytest

from mflco.model.params import TypicalSectionParameters
from mflco.aero.quasi_steady import QuasiSteady
from mflco.aero.qs_stall import QSStall, load_polar


def _section():
    """Any well-posed section: the equivalence identity is parameter-free."""
    return TypicalSectionParameters(a=-0.5, x_alpha=0.1, r_alpha_sq=0.5,
                                    omega_ratio=0.8, mu=150.0)


def test_linear_table_zero_trim_reduces_to_qs():
    """With Cl = 2*pi*alpha and trim 0, forces() must equal QuasiSteady to
    machine precision -- the machinery differs from QS by nothing but the
    table lookup. Tolerance 1e-14: observed 4e-17."""
    p = _section()
    lin_a = np.linspace(-60, 60, 241)
    qs = QuasiSteady(p, 0.0)
    st = QSStall(p, 0.0, lin_a, 2 * np.pi * np.radians(lin_a), trim_deg=0.0)
    rng = np.random.default_rng(0)
    for _ in range(100):
        y = rng.normal(0, 0.15, 4)
        U = rng.uniform(1.0, 5.0)
        assert st.forces(0, y, np.zeros(0), U) == pytest.approx(
            qs.forces(0, y, np.zeros(0), U), abs=1e-14)


def test_polar_peaks_at_the_rig_trim():
    """The bundled NACA 0020 / Re 1.6e5 polar peaks at 10-11 deg with
    Cl_max ~ 0.757 -- the rig's 10 deg trim sits AT static Cl_max. Tolerance
    0.02 on Cl: thickness interpolation of the Sheldahl-Klimas grid."""
    ag, cl = load_polar()
    # restrict to the pre-/near-stall lobe: the post-stall flat-plate recovery
    # climbs back to ~0.86 by 30 deg and would otherwise win the argmax
    i = (ag >= 0) & (ag <= 20)
    a_pk = ag[i][int(np.argmax(cl[i]))]
    assert 9.5 <= a_pk <= 11.5
    assert float(np.max(cl[i])) == pytest.approx(0.757, abs=0.02)


def test_slope_collapse_at_trim():
    """Local lift slope at 10 deg trim is ~0.63 /rad vs 2*pi -- the ~90%
    collapse that removes the QS flutter mechanism. Tolerance 0.05: central
    difference on the tabulated grid."""
    ag, cl = load_polar()
    st = QSStall(_section(), 0.0, ag, cl, trim_deg=10.0)
    assert st.slope_trim == pytest.approx(0.631, abs=0.05)
    assert st.slope_trim < 0.15 * 2 * np.pi


def test_linearisation_scales_with_trim_slope():
    """K_aero/C_aero are the QS matrices scaled by slope_trim/(2*pi) -- the
    linearisation about trim used by modal_analysis."""
    ag, cl = load_polar()
    p = _section()
    qs, st = QuasiSteady(p, 0.0), QSStall(p, 0.0, ag, cl, trim_deg=10.0)
    s = st.slope_trim / (2 * np.pi)
    assert st.K_aero(3.0) == pytest.approx(s * qs.K_aero(3.0), rel=1e-12)
    assert st.C_aero(3.0) == pytest.approx(s * qs.C_aero(3.0), rel=1e-12)


#verified:
# - QSStall(linear table, trim 0) == QuasiSteady to <1e-14 (observed 4e-17)
# - bundled polar peaks at 10-11 deg, Cl_max ~0.757: rig trims at static Cl_max
# - slope at trim ~0.63 /rad (90% collapse vs 2*pi) locked as a regression
# - trim-slope linearisation wiring for modal_analysis
#   (no-flutter / no-LCO probe result gated by examples/stage3_stall_probe.py)