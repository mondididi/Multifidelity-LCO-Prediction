"""Tests for analysis.py — eigenanalysis verification."""
import pytest

from mflco.model.analysis import coupled_natural_frequencies
from mflco.model.params import TypicalSectionParameters


def test_isogai_case_a_eigenfrequencies():
    """Wind-off coupled frequencies of Isogai Case A match Isogai (1979).
    Expected 0.713 and 5.34 in units of omega_alpha (Isogai 1979, reproduced in Yuan, Sandhu & Poirel, J. Aerospace Eng. 2021). 
    Tolerance 1% — Isogai's published values are 3 sig figs."""
    p = TypicalSectionParameters()   # defaults to Isogai Case A
    omegas, _ = coupled_natural_frequencies(p)
    
    assert omegas[0] == pytest.approx(0.713, rel=1e-2)
    assert omegas[1] == pytest.approx(5.34, rel=1e-2)


def test_uncoupled_limit_recovers_input_frequencies():
    """With x_alpha=0 the system decouples; frequencies are omega_h/omega_alpha and 1."""
    p = TypicalSectionParameters(x_alpha=0.0, omega_ratio=0.5, r_alpha_sq=1.0)
    omegas, _ = coupled_natural_frequencies(p)
    
    # With M, K both diagonal: K = diag(omega_ratio**2, r_alpha_sq), M = diag(1, r_alpha_sq).
    # Eigenvalues are diagonals of M^-1 K = (omega_ratio**2, 1). Square roots: (0.5, 1.0).
    assert omegas[0] == pytest.approx(0.5, abs=1e-12)
    assert omegas[1] == pytest.approx(1.0, abs=1e-12)

#verified:
# - Isogai Case A wind-off frequencies match published values within 1% (0.713 and 5.34)
# - In the uncoupled limit (x_alpha=0), natural frequencies match input parameters (omega_ratio and 1)