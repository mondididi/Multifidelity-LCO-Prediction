"""Tests for the Peters time-march path (Job B): the inflow dynamics aero_rhs,
the apparent-mass rework in eom.structural_rhs, and LCO behaviour.

The oracle is Job A <-> Job B agreement. The descriptor eigenproblem (Job A) and
the nonlinear time-march (Job B) are fully separate code paths -- a bug in one
will not surface in the other -- so a perturbed linear time-march must decay (or
grow) at the rate and frequency the least-damped eigenvalue predicts. If aero_rhs
forgot q'', used the bare M_s, or got a sign wrong, the marched rate would not
match the eigenvalue.
"""

import numpy as np
import pytest
from scipy.signal import find_peaks

from mflco.model.params import TypicalSectionParameters
from mflco.model.michigan_params import (
    calibrate_michigan, section_from_params, structural_zeta,
)
from mflco.aero.peters_finite import PetersFinite
from mflco.model.solver import integrate
from mflco.model.analysis import linearized_eigenvalues

cal = calibrate_michigan(zeta=structural_zeta())
p = section_from_params(cal)


def _least_damped_eig(P, U_star, N=6):
    """The structural eigenvalue closest to instability (the flutter mode)."""
    ev = linearized_eigenvalues(P, 0.0, U_star, PetersFinite(P, 0.0, N=N))
    osc = ev[ev.imag > 1e-8]
    zeta = -osc.real / np.abs(osc)
    return osc[np.argsort(zeta)][0]


def _fit_rate_and_freq(sig, tau):
    """Fit growth/decay rate sigma and damped frequency omega_d from a signal,
    using the |peak| envelope (log-linear slope) and the peak spacing."""
    pk, _ = find_peaks(np.abs(sig))
    pk = pk[np.abs(sig[pk]) > 1e-10]
    t_pk, a_pk = tau[pk], np.abs(sig[pk])
    sigma = np.polyfit(t_pk, np.log(a_pk), 1)[0]      # rate per tau
    omega_d = np.pi / np.mean(np.diff(t_pk))          # two peaks per period
    return sigma, omega_d


def test_solver_probe_shape():
    """aero_rhs returns (N,) so the solver's setup-time assertion holds. The
    degenerate stub used to return None (length-0 mismatch); this guards it."""
    aero = PetersFinite(p, M_inf=0.0, N=4)
    out = np.asarray(aero.aero_rhs(0.0, [0.0, 0.01, 0.0, 0.0], np.zeros(4), 2.0))
    assert out.shape == (4,)


@pytest.mark.parametrize("U_star", [2.0, 2.5])
def test_timemarch_decays_at_eigenvalue_rate(U_star):
    """Below flutter (~U*=3.67): a small perturbation decays at the least-damped
    eigenvalue's rate and frequency. This is the Job A <-> Job B cross-check."""
    lam = _least_damped_eig(p, U_star)
    sol = integrate(p, [0.0, 0.01, 0.0, 0.0], (0.0, 150.0),
                    aero=PetersFinite(p, 0.0, N=6), U_star=U_star,
                    t_eval=np.linspace(0.0, 150.0, 8000))
    sigma, omega_d = _fit_rate_and_freq(sol.y[1], sol.t)
    assert sigma < 0                                   # stable below flutter
    assert sigma == pytest.approx(lam.real, abs=6e-3)  # rate matches eigenvalue
    assert omega_d == pytest.approx(lam.imag, rel=0.03)


def test_timemarch_grows_above_flutter():
    """Above flutter the same perturbation grows, and the marched growth rate is
    positive and matches the (now unstable) eigenvalue's real part in sign."""
    U_star = 3.9
    lam = _least_damped_eig(p, U_star)
    assert lam.real > 0                                # eigenvalue says unstable
    sol = integrate(p, [0.0, 0.001, 0.0, 0.0], (0.0, 150.0),
                    aero=PetersFinite(p, 0.0, N=6), U_star=U_star,
                    t_eval=np.linspace(0.0, 150.0, 8000))
    sigma, _ = _fit_rate_and_freq(sol.y[1], sol.t)
    assert sigma > 0                                   # grows
    assert sigma == pytest.approx(lam.real, abs=6e-3)


def test_lco_is_bounded_with_cubic_spring():
    """With a hardening cubic pitch spring, a perturbation above the linear flutter
    speed saturates into a bounded limit cycle (not unbounded growth, not decay to
    zero). Confirms the nonlinear time-march closes the loop for LCO prediction."""
    p_lco = TypicalSectionParameters(
        a=cal.a, x_alpha=cal.x_alpha, r_alpha_sq=cal.r_alpha_sq,
        omega_ratio=cal.omega_ratio, mu=cal.mu, beta=40.0,
        zeta_h=cal.zeta, zeta_alpha=cal.zeta)
    sol = integrate(p_lco, [0.0, 0.02, 0.0, 0.0], (0.0, 600.0),
                    aero=PetersFinite(p_lco, 0.0, N=6), U_star=4.0,
                    t_eval=np.linspace(0.0, 600.0, 12000))
    early = sol.y[1][sol.t < 150]
    late = sol.y[1][sol.t > 450]
    amp_late = 0.5 * (late.max() - late.min())
    assert np.isfinite(late).all()                     # did not blow up
    assert amp_late > np.radians(0.5)                  # a real LCO, not decay
    assert amp_late < np.radians(30.0)                  # bounded (saturated)
    # amplitude has settled: late-window swing no longer growing materially
    assert late.max() <= early.max() * 5.0