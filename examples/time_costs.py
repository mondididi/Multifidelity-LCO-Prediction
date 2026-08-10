"""Stamp the per-query cost of every Python rung on THIS machine.

The cost metric of Sec. 3.5 is the wall time of one stability query at one
flight condition. The comparable primitive across rungs is a fixed-length
settle: integrate tau = 0 -> 650 from an 8 deg pitch kick at the amplitude
anchor speed (13.5 m/s), identical tolerances (LSODA, rtol 1e-7). Run this
once and put the three numbers in Table 1; the ONERA figure should agree
with the nominal characterisation record (probe_seconds / n_probes ~ 12 s).

Run:  PYTHONPATH=src python examples/time_costs.py     (~1 min)
"""
import time

import numpy as np
from scipy.integrate import solve_ivp

from mflco.model.params import TypicalSectionParameters
from mflco.model.michigan_params import calibrate_michigan, structural_zeta
from mflco.model.eom import structural_rhs
from mflco.aero.quasi_steady import QuasiSteady
from mflco.aero.peters_finite import PetersFinite
from mflco.aero.onera_stall import ONERAStall

U_MS, TAU_END, KICK = 13.5, 650.0, 8.0
cal = calibrate_michigan(zeta=structural_zeta())
p = TypicalSectionParameters(
    a=cal.a, x_alpha=cal.x_alpha, r_alpha_sq=cal.r_alpha_sq,
    omega_ratio=cal.omega_ratio, mu=cal.mu, beta=1.326, delta=0.0,
    zeta_h=cal.zeta, zeta_alpha=cal.zeta)

for name, aero in (("QS", QuasiSteady(p, 0.0)),
                   ("Peters N=6", PetersFinite(p, 0.0, N=6)),
                   ("ONERA (laws)", ONERAStall(p, 0.0, N=6))):
    n_a = getattr(aero, "n_aero_states", 0)
    y0 = np.zeros(4 + n_a)
    y0[1] = np.radians(KICK)
    t0 = time.perf_counter()
    solve_ivp(structural_rhs, (0.0, TAU_END), y0,
              args=(p, aero, float(cal.ms_to_ustar(U_MS))),
              method="LSODA", rtol=1e-7, atol=1e-9)
    print(f"  {name:14s} settle(650 tau) = {time.perf_counter()-t0:6.1f} s")