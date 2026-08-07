"""Stage 3 -- QS + static stall table: does static separation carry the LCO?

Structure held FIXED at the stage-2 build (delta = 0, beta = 1.326): the
ablation changes aerodynamics only, so any change in behaviour is
attributable to the stall table. Trim = 10 deg, the rig's mean AoA.

Findings this script demonstrates (first run 7/8/2026, reproduced in CI):
  1. The rig is trimmed AT the static Cl_max of its own airfoil: the NACA
     0020 polar at Re 1.6e5 peaks (~0.757) at 10-11 deg. Local lift slope at
     trim: ~0.63 /rad vs 2*pi = 6.28 -- a 90% collapse.
  2. Consequently the linearised QS+Stall section has NO flutter in range,
     and big-kick probes across the rig speed band find NO self-sustained
     LCO: static separation, quasi-steadily applied, produces nothing. The
     rung is priced out, and the split sharpens -- wake unsteadiness owns the
     flutter point (Peters, -0.31%), dynamic separation must own the fold
     (ONERA's rung).
Coherence note: a Peters + static-table hybrid is not attempted because the
inflow states assume attached 2*pi circulation dynamics; ONERA is the
consistent unsteadiness + separation combination.

Polar regeneration (provenance): clone github.com/SNL-WaterPower/CACTUS and
call build_polar_from_cactus(<repo>/test/Airfoil_Section_Data, out_csv) --
thickness interpolation of the Sheldahl-Klimas NACA 0018/0021 tables at
Re = 1.6e5. The bundled CSV in src/mflco/aero/data/ was built exactly so.

"""
import json
import os
import tempfile

import numpy as np

from mflco.model.params import TypicalSectionParameters
from mflco.model.michigan_params import calibrate_michigan, structural_zeta
from mflco.model.analysis import modal_analysis
from mflco.model.eom import structural_rhs
from mflco.model.fold_detector import classify_orbit
from mflco.aero.quasi_steady import QuasiSteady
from mflco.aero.qs_stall import QSStall, load_polar

TRIM_DEG   = 10.0      # rig mean AoA -- where the polar is sampled
BETA_S2    = 1.326     # stage-2 structural build, held fixed (delta = 0)
PROBE_U    = [8.0, 9.0, 10.0, 11.0, 11.85, 12.5, 13.19, 14.0]
PROBE_KICK = [8.0, -8.0, 15.0, -15.0]          # deg initial pitch
OUT_JSON   = os.path.join(tempfile.gettempdir(), "stage3_stall_probe.json")

cal = calibrate_michigan(zeta=structural_zeta())


def _section():
    return TypicalSectionParameters(
        a=cal.a, x_alpha=cal.x_alpha, r_alpha_sq=cal.r_alpha_sq,
        omega_ratio=cal.omega_ratio, mu=cal.mu, beta=BETA_S2, delta=0.0,
        zeta_h=cal.zeta, zeta_alpha=cal.zeta)


def build_polar_from_cactus(cactus_data_dir, out_csv, re_tag="1.6e5"):
    """Regenerate the NACA 0020 polar from the CACTUS Sheldahl-Klimas tables."""
    def _block(path):
        rows, on = [], False
        for line in open(path):
            if line.startswith("Reynolds Number:"):
                on = line.split(":")[1].strip() == re_tag
                continue
            if on:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        rows.append((float(parts[0]), float(parts[1])))
                    except ValueError:
                        pass
                elif rows:
                    break
        return np.array(rows)
    t18 = _block(os.path.join(cactus_data_dir, "NACA_0018.dat"))
    t21 = _block(os.path.join(cactus_data_dir, "NACA_0021.dat"))
    grid = np.union1d(t18[:, 0], t21[:, 0])
    cl20 = (np.interp(grid, t18[:, 0], t18[:, 1]) / 3.0
            + 2.0 * np.interp(grid, t21[:, 0], t21[:, 1]) / 3.0)
    with open(out_csv, "w") as f:
        f.write("# NACA 0020, Re=1.6e5. Thickness interpolation of "
                "Sheldahl-Klimas\n# (SAND80-2114) NACA 0018/0021 tables via "
                "Sandia CACTUS data files.\n# alpha_deg, Cl\n")
        for a, c in zip(grid, cl20):
            f.write(f"{a:.2f},{c:.6f}\n")


def _flutter_or_none(aero, us_lo=0.3, us_hi=6.0, n=400):
    prev, Uf = None, None
    for u in np.linspace(us_lo, us_hi, n):
        try:
            _, d = modal_analysis(aero.params, 0.0, u, aero)
        except ValueError:
            prev = None
            continue
        m = float(np.min(d))
        if prev is not None and prev > 0.0 >= m:
            Uf = u
            break
        prev = m
    return None if Uf is None else float(cal.ustar_to_ms(Uf))


# --- 1. machinery check: linear table + zero trim must equal QS -------------------
p = _section()
lin_a = np.linspace(-60, 60, 241)
qs = QuasiSteady(p, 0.0)
st_lin = QSStall(p, 0.0, lin_a, 2 * np.pi * np.radians(lin_a), trim_deg=0.0)
rng = np.random.default_rng(0)
err = max(float(np.max(np.abs(qs.forces(0, y, np.zeros(0), U)
                              - st_lin.forces(0, y, np.zeros(0), U))))
          for y, U in ((rng.normal(0, 0.15, 4), rng.uniform(1, 5))
                       for _ in range(200)))
print(f"equivalence |QS - QSStall(linear table, trim 0)|_max = {err:.3e}")

# --- 2. the real polar at trim ----------------------------------------------------
ag, cl = load_polar()
aero = QSStall(p, 0.0, ag, cl, trim_deg=TRIM_DEG)
i_max = int(np.argmax(cl[(np.abs(ag) <= 30)]))
print(f"polar at trim {TRIM_DEG} deg: Cl = {aero._cl_trim:.4f}, "
      f"slope = {aero.slope_trim:.3f} /rad  (2*pi = {2*np.pi:.3f})")

# --- 3. linear flutter attempt ----------------------------------------------------
Uf_qs = _flutter_or_none(QuasiSteady(p, 0.0))
Uf_st = _flutter_or_none(aero)
print(f"linear flutter, pure QS      : "
      f"{'none' if Uf_qs is None else f'{Uf_qs:.3f} m/s'}")
print(f"linear flutter, QS+Stall({TRIM_DEG:.0f}) : "
      f"{'none' if Uf_st is None else f'{Uf_st:.3f} m/s'}")

# --- 4. finite-amplitude LCO probe ------------------------------------------------
print("finite-amplitude LCO probe (kick -> sustained LCO or decay):")
found = []
for U in PROBE_U:
    hits = []
    for a0 in PROBE_KICK:
        y0 = np.zeros(4)
        y0[1] = np.radians(a0)
        amp, _, k, fl = classify_orbit(
            structural_rhs, (p, aero, float(cal.ms_to_ustar(U))), y0)
        if amp is not None:
            hits.append((a0, round(amp, 2)))
    tag = " / ".join(f"kick{a0:+.0f} -> LCO {am:.1f} deg" for a0, am in hits) \
        or "all kicks DEAD"
    print(f"  U={U:6.2f}  {tag}", flush=True)
    if hits:
        found.append((U, hits))

result = dict(model="QS+Stall", trim=TRIM_DEG, beta=BETA_S2,
              equivalence_err=err,
              slope_trim=round(aero.slope_trim, 4),
              Uf_qs=Uf_qs, Uf_stall=Uf_st,
              lco_found=found or "none",
              probes=dict(U=PROBE_U, kick_deg=PROBE_KICK))
json.dump(result, open(OUT_JSON, "w"), indent=1)
print(f"\nresult: {'LCOs found -- investigate' if found else 'no self-sustained LCO at any probed speed/kick: static stall alone does not produce the oscillation'}")
print(f"json -> {OUT_JSON}")