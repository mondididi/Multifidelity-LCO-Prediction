"""Barton (AeroelasticCBC.jl src/model.jl) parameters in mflco form.

Second-rig parameter module, sibling of michigan_params.py. Purpose: verify
the structural machinery with ZERO free parameters. Every number below is
published in Barton's model.jl (Tartaruga/Bristol rig, NACA 0015, chord
0.30 m); the mapping to the TypicalSectionParameters nondimensional groups is
closed-form; wind-off coupled frequencies can then be checked against the raw
DIMENSIONAL eigenproblem computed independently in this module. Agreement to
machine precision proves the parameter container, M/K assembly and scaling
carry no hidden assumptions. See examples/stage0_barton_test.py.

Mapping (heave carries the TOTAL mass m_T; only the wing m_w pitches):
    x_alpha_eff = m_w * x_alpha / m_T
    r_alpha_sq  = I_alpha / (m_T b^2)
    omega_alpha = sqrt(k_alpha / I_alpha);  omega_h = sqrt(k_h / m_T)
    omega_ratio = omega_h / omega_alpha
    mu          = m_T / (pi rho b^2)
    beta        = k_alpha3 / k_alpha            (nondimensional cubic)
    delta       = k_alpha2 / (3 k_alpha3)       (Barton's FREE quadratic is
                                                 exactly the constrained
                                                 trim-offset form of params.py)
    zeta_h      = c_h / (2 m_T omega_h)
    zeta_alpha  = c_alpha / (2 I_alpha omega_alpha)
    U [m/s]     = U* * b * omega_alpha

Known answers for case 2 (the published experimental campaign):
    full-Sears standalone port : Hopf 24.0 m/s, fold 16.81 m/s, SUBCRITICAL
    Tartaruga rig (IFASD 2019) : Hopf ~24 m/s, fold ~16 m/s, subcritical
"""

import numpy as np

from mflco.model.params import TypicalSectionParameters

# --- Barton's published numbers (model.jl, verbatim) -----------------------------
BASE = dict(b=0.15, a=-0.5, rho=1.204, m_w=5.3, m_T=16.9)
CASES = {
    1: dict(I_a=0.1724, c_a=0.5628, c_h=14.5756, k_a=54.1162,
            k_a2=751.6, k_a3=5006.7, k_h=3529.4, x_a=0.24),
    2: dict(I_a=0.1726, c_a=1.0338, c_h=15.4430, k_a=60.291,
            k_a2=774.7, k_a3=3490.7, k_h=3317.8, x_a=0.234),
}
KNOWN_CASE2 = dict(hopf_ms=24.0, fold_ms=16.81,       # full-Sears port
                   exp_hopf_ms=24.0, exp_fold_ms=16.0)  # Tartaruga rig, approx


class BartonCal:
    """Closed-form nondimensionalisation of one Barton case."""

    def __init__(self, case):
        P = dict(BASE, **CASES[case])
        self.case = case
        self.b, self.a = P["b"], P["a"]
        self.omega_alpha = np.sqrt(P["k_a"] / P["I_a"])
        self.omega_h = np.sqrt(P["k_h"] / P["m_T"])
        self.omega_ratio = self.omega_h / self.omega_alpha
        self.x_alpha = P["m_w"] * P["x_a"] / P["m_T"]
        self.r_alpha_sq = P["I_a"] / (P["m_T"] * P["b"] ** 2)
        self.mu = P["m_T"] / (np.pi * P["rho"] * P["b"] ** 2)
        self.beta = P["k_a3"] / P["k_a"]
        self.delta = P["k_a2"] / (3.0 * P["k_a3"])
        self.zeta_h = P["c_h"] / (2.0 * P["m_T"] * self.omega_h)
        self.zeta_alpha = P["c_a"] / (2.0 * P["I_a"] * self.omega_alpha)
        self._P = P

    def ms_to_ustar(self, U):
        return U / (self.b * self.omega_alpha)

    def ustar_to_ms(self, Us):
        return Us * self.b * self.omega_alpha

    def section(self):
        return TypicalSectionParameters(
            a=self.a, x_alpha=self.x_alpha, r_alpha_sq=self.r_alpha_sq,
            omega_ratio=self.omega_ratio, mu=self.mu,
            beta=self.beta, delta=self.delta,
            zeta_h=self.zeta_h, zeta_alpha=self.zeta_alpha)

    def windoff_dimensional_hz(self):
        """Coupled wind-off frequencies from the RAW dimensional matrices --
        computed with no shared mflco code: the zero-assumption reference."""
        P = self._P
        S = P["m_w"] * P["x_a"] * P["b"]
        M = np.array([[P["m_T"], S], [S, P["I_a"]]])
        K = np.array([[P["k_h"], 0.0], [0.0, P["k_a"]]])
        lam = np.sort(np.linalg.eigvals(np.linalg.solve(M, K)).real)
        return np.sqrt(lam) / (2.0 * np.pi)