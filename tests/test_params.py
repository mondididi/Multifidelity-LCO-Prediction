"""Tests for TypicalSectionParameters.""" #file must start with test_ to be discovered by pytest
import numpy as np
import pytest

from mflco.model.params import TypicalSectionParameters


def test_default_params_match_isogai_case_a():
    """Default constructor returns Isogai Case A values."""
    p = TypicalSectionParameters()
    assert p.a           == -2.0
    assert p.x_alpha     == 1.8
    assert p.r_alpha_sq  == 3.48
    assert p.omega_ratio == 1.0
    assert p.mu          == 60.0
    assert p.beta        == 0.0
    assert p.zeta_h      == 0.0
    assert p.zeta_alpha  == 0.0


def test_mass_matrix_isogai_values():
    """Mass matrix equals Isogai Case A values."""
    p = TypicalSectionParameters()
    expected = np.array([[1.0, 1.8],
                         [1.8, 3.48]])
    np.testing.assert_array_equal(p.mass_matrix(), expected)


def test_mass_matrix_symmetric():
    """Mass matrix off-diagonals are equal (symmetric)."""
    p = TypicalSectionParameters()
    M = p.mass_matrix()
    assert M[0, 1] == M[1, 0]


def test_damping_matrix_zero_when_undamped():
    """Damping matrix is all zeros when zeta_h = zeta_alpha = 0 (default)."""
    p = TypicalSectionParameters()
    C = p.damping_matrix()
    assert np.all(C == 0)


def test_stiffness_matrix_linear_when_beta_zero():
    """With beta = 0, K does not depend on alpha."""
    p = TypicalSectionParameters()   # beta defaults to 0
    np.testing.assert_array_equal(p.stiffness_matrix(0.0),
                                  p.stiffness_matrix(0.3))


def test_stiffness_matrix_hardens_when_beta_positive():
    """With beta = 3, K[1,1] grows with alpha² per (1 + beta alpha²)."""
    p = TypicalSectionParameters(beta=3.0)
    K = p.stiffness_matrix(0.3)
    # Expected: r_alpha² · (1 + beta alpha²) = 3.48 · (1 + 3·0.09) = 4.4196
    assert K[1, 1] == pytest.approx(4.4196)


#verified:
# - mass matrix matches Isogai Case A values
# - mass matrix is symmetric
# - damping matrix is zero when undamped
# - stiffness matrix is linear when beta = 0
# - stiffness matrix hardens with alpha when beta > 0