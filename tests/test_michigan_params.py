"""
test_michigan_params.py
=======================
Gate the structural calibration: the calibrated section must reproduce the
measured U=0 modes (5.3 / 6.2 Hz). Flips green once michigan_params is correct.

The `xfail` markers document the contract; remove `strict`/the marker once
michigan_section() is wired to the real mflco TypicalSection and you want the
green check to be load-bearing in CI.
"""

import numpy as np
import pytest
from mflco.model.michigan_params import michigan_section
from mflco.model.michigan_params import (
    calibrate_michigan, F_PLUNGE_HZ, F_PITCH_HZ,
    M_KG, S_ALPHA_KGM, B_M,
)

TARGET = (F_PLUNGE_HZ, F_PITCH_HZ)


def _coupled_modes_hz(p):
    """U=0 coupled modes from the *effective dimensional* structure (Hz)."""
    M = np.array([[M_KG, S_ALPHA_KGM], [S_ALPHA_KGM, p.Ia_eff]])
    K = np.array([[p.Kh_eff, 0.0], [0.0, p.Ka_eff]])
    w2 = np.sort(np.linalg.eigvals(np.linalg.solve(M, K)).real)
    return np.sqrt(w2) / (2 * np.pi)


def test_u0_modes_match_measured():
    """Calibrated section reproduces 5.3 / 6.2 Hz at U=0."""
    p = calibrate_michigan(omega_ratio=1.0)
    f_lo, f_hi = _coupled_modes_hz(p)
    assert f_lo == pytest.approx(TARGET[0], abs=1e-3)
    assert f_hi == pytest.approx(TARGET[1], abs=1e-3)


def test_mode_ratio_matches():
    """Catches a bad Fig-6 reading: the *ratio* must be reproduced too."""
    p = calibrate_michigan(omega_ratio=1.0)
    f_lo, f_hi = _coupled_modes_hz(p)
    assert (f_hi / f_lo) == pytest.approx(TARGET[1] / TARGET[0], rel=1e-4)


def test_geometry_held_fixed():
    """a, x_alpha, mu must come from geometry, not be tuned."""
    p = calibrate_michigan()
    assert p.a == pytest.approx(-0.5)
    assert p.x_alpha == pytest.approx(0.11875, rel=1e-4)
    assert p.mu == pytest.approx(141.7, rel=2e-3)


def test_velocity_scale_robust():
    """U-scale (b*omega_alpha) is ~3.5-3.6 m/s/U* across plausible readings."""
    for f_lo, f_hi in [(5.3, 6.2), (5.0, 6.5), (4.8, 6.8)]:
        p = calibrate_michigan(f_lo, f_hi, omega_ratio=1.0)
        assert 3.3 < p.b * p.omega_alpha < 3.8


@pytest.mark.parametrize("wr", [0.95, 1.0, 1.05])
def test_family_members_all_reproduce_modes(wr):
    """Any in-family omega_ratio pin still hits the measured modes."""
    p = calibrate_michigan(omega_ratio=wr)
    f_lo, f_hi = _coupled_modes_hz(p)
    assert (f_lo, f_hi) == pytest.approx(TARGET, abs=1e-3)


@pytest.mark.xfail(reason="wire michigan_section() to mflco.TypicalSection first",
                   strict=False)
def test_section_object_built():
    """Once wired, michigan_section() returns a real section, not the params."""

    sec = michigan_section()
    # e.g. assert sec.linearized_eigenvalues(0.0) gives U=0 modes 5.3/6.2 Hz
    assert hasattr(sec, "system_matrix")