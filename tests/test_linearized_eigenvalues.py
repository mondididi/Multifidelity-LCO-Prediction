"""
Tests for linearized eigenvalue analysis (4×4 first-order method).

Verifies:
- Structural correctness of the system matrix Q. (4x4)
- Agreement with the 2×2 undamped special case (wind-off, ζ=0).
    - At U*=0, ζ=0, no aero: all eigenvalues should be pure imaginary (no damping).
    - The 4×4 method's |Im λ| should match the 2×2 method's frequencies.
- Conjugate-pair structure of eigenvalues.
- Nonlinear stiffness effect: frequencies shift with α_eq when β > 0.
"""

import numpy as np
import pytest
from mflco.model.params import TypicalSectionParameters
from mflco.model.analysis import (
    system_matrix,
    linearized_eigenvalues,
    undamped_natural_frequencies,
)


def test_system_matrix_shape_and_structure():
    """Q is 4×4 with the expected block structure: zeros top-left, identity top-right."""
    p = TypicalSectionParameters()
    Q = system_matrix(p)
    
    # Shape
    assert Q.shape == (4, 4), f"Expected Q shape (4,4), got {Q.shape}"
    
    # Top-left 2×2 block should be all zeros
    assert np.allclose(Q[:2, :2], np.zeros((2, 2))), \
        "Top-left block of Q should be zero (positions don't depend directly on themselves)"
    
    # Top-right 2×2 block should be identity (kinematic identity: dq/dt = q_dot)
    assert np.allclose(Q[:2, 2:], np.eye(2)), \
        "Top-right block of Q should be identity (the trivial kinematic identity)"


def test_eigenvalues_pure_imaginary_at_wind_off():
    """At U*=0, ζ=0, no aero: all eigenvalues should have Re(λ) ≈ 0 (no damping)."""
    p = TypicalSectionParameters()
    eigs = linearized_eigenvalues(p)
    
    # All real parts should be near zero (within numerical noise)
    max_real = np.max(np.abs(eigs.real))
    assert max_real < 1e-10, \
        f"Wind-off eigenvalues should be pure imaginary; max |Re(λ)| = {max_real:.2e}"


def test_eigenvalues_match_undamped_special_case():
    """At wind-off, the 4×4 method's |Im λ| should match the 2×2 method's frequencies."""
    p = TypicalSectionParameters()
    
    # New 4×4 method
    eigs = linearized_eigenvalues(p)
    # Take only the positive-imaginary eigenvalues (one of each conjugate pair)
    positive_freqs = sorted(np.abs(eigs.imag[eigs.imag > 0]))
    
    # Old 2×2 method
    omegas, _ = undamped_natural_frequencies(p)
    expected_freqs = sorted(omegas)
    
    assert np.allclose(positive_freqs, expected_freqs, rtol=1e-10), \
        f"4×4 method gave {positive_freqs}, 2×2 method gave {expected_freqs}"


def test_eigenvalues_come_in_conjugate_pairs():
    """For a real Q matrix, complex eigenvalues must appear as conjugate pairs."""
    p = TypicalSectionParameters()
    eigs = linearized_eigenvalues(p)
    
    # For each eigenvalue, its complex conjugate must also appear in the list
    # (within numerical tolerance)
    for e in eigs:
        conjugate_present = any(np.isclose(e.conjugate(), other, atol=1e-10) for other in eigs)
        assert conjugate_present, \
            f"Eigenvalue {e} has no conjugate partner in {eigs}"


def test_frequencies_shift_with_alpha_eq_when_beta_positive():
    """
    With cubic spring β > 0, linearizing K around α_eq > 0 stiffens the pitch DOF
    (k_α · (1 + β·α_eq²) > k_α), so frequencies should shift.
    
    confirms the nonlinear stiffness is wired correctly into the eigenvalue analysis, 
    not just into time-history integration.
    """
    # Use β > 0 so the cubic spring is active
    p = TypicalSectionParameters(beta=3.0)
    
    # Frequencies at α_eq = 0 (linear stiffness)
    eigs_linear = linearized_eigenvalues(p, alpha_eq=0.0)
    freqs_linear = sorted(np.abs(eigs_linear.imag[eigs_linear.imag > 0]))
    
    # Frequencies at α_eq = 0.2 rad (cubic stiffness active)
    eigs_nonlinear = linearized_eigenvalues(p, alpha_eq=0.2)
    freqs_nonlinear = sorted(np.abs(eigs_nonlinear.imag[eigs_nonlinear.imag > 0]))
    
    # The pitch-dominated mode (higher frequency) should stiffen most
    # because β affects K[1,1] specifically
    assert freqs_nonlinear[1] > freqs_linear[1], \
        f"With β>0, frequency at α_eq=0.2 ({freqs_nonlinear[1]:.4f}) should exceed " \
        f"frequency at α_eq=0 ({freqs_linear[1]:.4f})"