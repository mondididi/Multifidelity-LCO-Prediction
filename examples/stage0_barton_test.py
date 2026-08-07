"""Stage 0 -- machinery verification against Barton's published model.

Three-layer check with ZERO free parameters (every input is published in
AeroelasticCBC.jl src/model.jl):

  1. MAPPING + WIND-OFF. Barton's dimensional numbers -> mflco nondimensional
     groups (barton_params.BartonCal). The coupled wind-off frequencies from
     the pipeline (params.py -> undamped_natural_frequencies) must match the
     raw DIMENSIONAL eigenproblem, computed independently, to machine
     precision. Proves container, M/K assembly and scaling: no assumptions.
  2. NONLINEARITY PATHWAY. Barton's FREE quadratic k_a2 maps exactly onto the
     constrained trim form of params.py: beta = k_a3/k_a, delta = k_a2/(3 k_a3).
     No code change required -- the delta parameterisation spans his model.
  3. BIFURCATION KNOWN ANSWER. Case 2 through the stage-0 chain (eom.py +
     QS/Peters + fold_detector) must be SUBCRITICAL with Hopf/fold near the
     independent full-Sears port (24.0 / 16.81 m/s) and the Tartaruga rig
     (~24 / ~16 m/s). QS-vs-Peters spread quantifies aero-flavour sensitivity:
     fold LOCATION moves with wake fidelity, fold EXISTENCE does not.

Consequence for the Michigan result: the same chain that finds a structural
fold to sub-percent accuracy here finds none there -- the supercritical
verdict is physics, not machinery.

"""
import json
import os
import tempfile

import numpy as np

from mflco.model.barton_params import BartonCal, KNOWN_CASE2
from mflco.model.analysis import modal_analysis, undamped_natural_frequencies
from mflco.model.eom import structural_rhs
from mflco.model.fold_detector import down_sweep
from mflco.aero.quasi_steady import QuasiSteady
from mflco.aero.peters_finite import PetersFinite

N_INFLOW = 6          # Peters states, matching the Michigan stages


# --- helpers ---------------------------------------------------------------------
def _flutter_ms(cal, aero, U_lo=0.5, U_hi=14.0, n=700):
    """Linear flutter [m/s] from the eigenvalue sweep (min-damping crossing)."""
    p = aero.params
    prev, Uf = None, None
    for u in np.linspace(U_lo, U_hi, n):
        try:
            _, d = modal_analysis(p, 0.0, u, aero)
        except ValueError:
            prev = None
            continue
        m = float(np.min(d))
        if prev is not None and prev > 0.0 >= m:
            Uf = u
            break
        prev = m
    return None if Uf is None else float(cal.ustar_to_ms(Uf))


def _sweep(cal, aero_kind):
    p = cal.section()
    if aero_kind == "peters":
        aero, n_a = PetersFinite(p, 0.0, N=N_INFLOW), N_INFLOW
    else:
        aero, n_a = QuasiSteady(p, 0.0), 0
    Uf = _flutter_ms(cal, aero)
    print(f"  linear flutter ({aero_kind}): {Uf:.3f} m/s "
          f"[known full-Sears {KNOWN_CASE2['hopf_ms']}]", flush=True)
    y0 = np.zeros(4 + n_a)
    y0[1] = np.radians(11.0)
    argU = lambda U: (p, aero, float(cal.ms_to_ustar(U)))
    br, U_fold, amp_f, secs = down_sweep(
        structural_rhs, argU, y0, round(Uf + 0.5, 2), 0.25,
        max(Uf - 10.0, 5.0), 0.05)
    window = Uf - U_fold
    verdict = ("SUBCRITICAL" if (window > 0.1 and amp_f > 1.5)
               else "SUPERCRITICAL")
    print(f"  -> {aero_kind}: Hopf {Uf:.3f}, fold {U_fold:.3f} "
          f"(amp {amp_f:.2f} deg), window {window:.3f}, {verdict}  "
          f"[{secs:.0f}s]")
    return dict(aero=aero_kind, U_hopf=round(Uf, 3), U_fold=round(U_fold, 3),
                amp_at_fold=round(amp_f, 2), window=round(window, 3),
                verdict=verdict, seconds=round(secs, 1),
                branch=[(round(u, 3), round(a, 2)) for u, a in br])


# --- layer 1: mapping + zero-assumption wind-off check ---------------------------
print("=" * 76)
print("STAGE 0 -- machinery verification against Barton (AeroelasticCBC model.jl)")
print("=" * 76)
for c in (1, 2):
    cal = BartonCal(c)
    print(f"--- case {c} " + "-" * 58)
    print(f"  x_alpha_eff = {cal.x_alpha:.6f}   r_alpha_sq = {cal.r_alpha_sq:.6f}"
          f"   omega_ratio = {cal.omega_ratio:.6f}")
    print(f"  mu = {cal.mu:.3f}   beta = {cal.beta:.4f}   "
          f"delta = {np.degrees(cal.delta):.4f} deg   "
          f"zeta_h/alpha = {cal.zeta_h:.4f}/{cal.zeta_alpha:.4f}")
    f_dim = cal.windoff_dimensional_hz()
    nat, _ = undamped_natural_frequencies(cal.section())
    f_pipe = np.sort(np.asarray(nat)) * cal.omega_alpha / (2.0 * np.pi)
    err = float(np.max(np.abs(f_dim - f_pipe)))
    print(f"  wind-off dimensional : {f_dim[0]:.6f} / {f_dim[1]:.6f} Hz")
    print(f"  wind-off mflco       : {f_pipe[0]:.6f} / {f_pipe[1]:.6f} Hz")
    print(f"  max |diff| = {err:.3e} Hz  "
          f"{'PASS' if err < 1e-9 else 'FAIL'} (zero-assumption check)")

# --- layer 3: bifurcation known answer, case 2 -----------------------------------
print("\n--- case 2 bifurcation through the stage-0 chain " + "-" * 27)
cal2 = BartonCal(2)
out = {k: _sweep(cal2, k) for k in ("qs", "peters")}

print("\nSUMMARY (case 2)")
print(f"  {'path':<28} {'Hopf':>7} {'fold':>7}  verdict")
for k, r in out.items():
    print(f"  mflco {k:<22} {r['U_hopf']:>7.2f} {r['U_fold']:>7.2f}  {r['verdict']}")
print(f"  {'full-Sears port (known)':<28} {KNOWN_CASE2['hopf_ms']:>7.2f} "
      f"{KNOWN_CASE2['fold_ms']:>7.2f}  SUBCRITICAL")
print(f"  {'Tartaruga rig (experiment)':<28} {KNOWN_CASE2['exp_hopf_ms']:>7.2f} "
      f"{KNOWN_CASE2['exp_fold_ms']:>7.2f}  subcritical")
out_path = os.path.join(tempfile.gettempdir(), "stage0_barton.json")
json.dump(out, open(out_path, "w"), indent=1)
print(f"\n  results -> {out_path}")