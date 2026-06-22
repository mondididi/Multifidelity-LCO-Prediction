"""Stage 2 (Peters finite-state inflow) VGBF on the Michigan section.

Same calibrated structure as the Stage 1 script (U=0 modes pinned to 5.3/6.2 Hz),
but the quasi-steady aero is replaced by Peters' N-state inflow -- the ONLY change
from stage1_michigan_eigenvalue_sweep.py is the aero object (QuasiSteady -> Peters)
plus the N knob. Everything downstream (descriptor eigenproblem, mode reduction,
flutter detection) is dispatched automatically because PetersFinite reports
n_aero_states > 0, so modal_analysis takes the (4+N) generalized path.

Result (vs experiment 13.19 m/s, fold 11.85):
    Stage 1 quasi-steady              ~ 9.6  m/s   (-27%, missing wake lag C(k))
    Stage 2 Peters, full (N>=5)       ~13.1  m/s   (-1%, wake lag + apparent mass)

The wake lag (Peters C(k)) recovers most of the QS deficit; the apparent
mass/damping closes the rest. Both must be present and correctly signed.
"""
import numpy as np
import matplotlib.pyplot as plt

from mflco.model.michigan_params import (
    calibrate_michigan, section_from_params, structural_zeta,
)
from mflco.aero.peters_finite import PetersFinite           # <-- the swap
from mflco.model.analysis import modal_analysis, undamped_natural_frequencies

# experimental landmarks (Garcia Perez et al., dimensional)
F_U0_HZ      = (5.3, 6.2)
U_FLUTTER_MS = 13.19
U_FOLD_MS    = 11.85

N_INFLOW     = 6               # Peters inflow states (flutter is converged by N~5)
ZETA_STRUCT  = structural_zeta()

cal = calibrate_michigan(zeta=ZETA_STRUCT)
p   = section_from_params(cal)
aero = PetersFinite(p, M_inf=0.0, N=N_INFLOW)               # incompressible rig

u_sweep_val = np.linspace(0.0, 4.5, 200)                    # non-dim U*

freq = np.full((len(u_sweep_val), 2), np.nan)
damp = np.full((len(u_sweep_val), 2), np.nan)
for i, U_star in enumerate(u_sweep_val):
    try:
        freq[i], damp[i] = modal_analysis(p, alpha_eq=0.0, U_star=U_star, aero=aero)
    except ValueError:
        continue

f_hz = freq * cal.omega_alpha / (2 * np.pi)
u_ms = cal.ustar_to_ms(u_sweep_val)

# U*=0 oracle: turn apparent mass off so the wind-off modes are the pure
# structural 5.3/6.2 Hz (apparent mass shifts them ~0.5% and would fail allclose).
aero_oracle = PetersFinite(p, M_inf=0.0, N=N_INFLOW, include_apparent_mass=False)
f0, _ = modal_analysis(p, 0.0, 0.0, aero_oracle)
nat, _ = undamped_natural_frequencies(p)
assert np.allclose(np.sort(f0), np.sort(nat), rtol=1e-3), (f0, nat)

# flutter = lowest U* where a structural mode's damping crosses + -> -
min_damp = np.nanmin(damp, axis=1)
U_flutter = None
for i in range(1, len(u_sweep_val)):
    a, b = min_damp[i - 1], min_damp[i]
    if np.isfinite(a) and np.isfinite(b) and a > 0.0 >= b:
        x0, x1 = u_sweep_val[i - 1], u_sweep_val[i]
        U_flutter = x0 - a * (x1 - x0) / (b - a)
        break
U_flutter_ms = None if U_flutter is None else float(cal.ustar_to_ms(U_flutter))

fig, (ax_f, ax_d) = plt.subplots(2, 1, figsize=(6, 7))
ax_f.plot(u_ms, f_hz[:, 0], label=f"Peters N={N_INFLOW} mode 1")
ax_f.plot(u_ms, f_hz[:, 1], label=f"Peters N={N_INFLOW} mode 2")
ax_f.plot([0, 0], F_U0_HZ, "ks", ms=6, label="exp. U=0 (Fig. 6)")
ax_f.set_ylabel("frequency [Hz]")

ax_d.axhline(0.0, color="k", lw=0.8)
ax_d.plot(u_ms, damp[:, 0]); ax_d.plot(u_ms, damp[:, 1])
ax_d.set_ylabel(r"damping ratio $\zeta$"); ax_d.set_xlabel("airspeed $U$ [m/s]")

for ax in (ax_f, ax_d):
    ax.axvline(U_FLUTTER_MS, color="r", ls="--", lw=1.0, label=f"exp. flutter {U_FLUTTER_MS} m/s")
    ax.axvline(U_FOLD_MS,    color="r", ls=":",  lw=1.0, label=f"exp. fold {U_FOLD_MS} m/s")

if U_flutter_ms is not None:
    for ax in (ax_f, ax_d):
        ax.axvline(U_flutter_ms, color="C2", ls="-.", lw=1.0,
                   label=f"Peters flutter {U_flutter_ms:.1f} m/s")
    ax_f.set_title(f"Stage 2 Peters flutter at U = {U_flutter_ms:.1f} m/s  (exp. {U_FLUTTER_MS})")
else:
    ax_f.set_title("Stage 2 Peters: no flutter crossing in swept range")

ax_f.legend(fontsize=8)
plt.tight_layout()
print(f"U*=0 oracle (apparent mass off): {np.sort(f0) * cal.omega_alpha / (2*np.pi)} Hz  (target 5.3 / 6.2)")
print(f"Stage 2 Peters N={N_INFLOW} flutter: "
      f"{U_flutter_ms:.2f} m/s" if U_flutter_ms else "no crossing"
      f"  (exp. {U_FLUTTER_MS}, fold {U_FOLD_MS})")

plt.show()