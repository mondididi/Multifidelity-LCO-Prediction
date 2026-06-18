"""Stage 1 (quasi-steady) VGBF on the Michigan section — experimental validation.

Calibrated Michigan structure (U=0 modes pinned to 5.3 / 6.2 Hz) + quasi-steady
aero. Sweep U*, build the system matrix, take eigenvalues, plot frequency &
damping vs *dimensional* speed (Hz, m/s) so it overlays García Pérez Fig. 6.
"""
# NOTE: with zeta=0 the pitch mode is unstable from U*~=0 (no structural damping
# + no C(k) lag = QS over-destabilises), so "no flutter crossing" is correct,
# not a bug. Back out zeta from the recovery rates for a finite QS flutter speed.

from mflco.model.michigan_params import calibrate_michigan, section_from_params
from mflco.aero.quasi_steady import QuasiSteady
from mflco.model.analysis import modal_analysis, undamped_natural_frequencies
import numpy as np
import matplotlib.pyplot as plt

# experimental landmarks (García Pérez, dimensional)
F_U0_HZ      = (5.3, 6.2)    # Fig. 6 U=0 intercepts
U_FLUTTER_MS = 13.19         # linear flutter onset
U_FOLD_MS    = 11.85         # subcritical fold

# initialise — Michigan section (calibrated); keep cal for the Hz / m/s scales
cal = calibrate_michigan()                   # omega_ratio=1.0, zeta=0.0 defaults, output as container
p   = section_from_params(cal)               # calibrated section (beta=0, linear), output as section for aeromodel, (arg: container)
qs  = QuasiSteady(p, M_inf=0.0)              # incompressible; swap to Peters later

u_sweep_val = np.linspace(0.0, 4.5, 90)      # non-dim U* (flutter ~3.7; 4.5 brackets it)

freq = np.full((len(u_sweep_val), 2), np.nan)   # (n_speeds, 2 modes)
damp = np.full((len(u_sweep_val), 2), np.nan)

# steps 1-7: sweep, eigenvalues, reduce to (freq, zeta), store, skip coalescence
for i, U_star in enumerate(u_sweep_val):
    try:
        freq[i], damp[i] = modal_analysis(p, alpha_eq=0.0, U_star=U_star, aero=qs)
    except ValueError:        # past coalescence -> non-oscillatory; leave NaN
        continue

# convert to dimensional units for the Fig. 6 overlay
f_hz = freq * cal.omega_alpha / (2 * np.pi)  # nondim omega/omega_alpha -> Hz
u_ms = cal.ustar_to_ms(u_sweep_val)          # U* -> m/s

# step 10: U*=0 oracle (aero off at zero speed -> must equal structural freqs)
nat, _ = undamped_natural_frequencies(p)
assert np.allclose(np.sort(freq[0]), np.sort(nat)), (freq[0], nat)

# step 8: flutter = lowest U* where any mode's damping crosses + -> -
min_damp = np.nanmin(damp, axis=1)
U_flutter = None
for i in range(1, len(u_sweep_val)):
    a, b = min_damp[i-1], min_damp[i]
    if np.isfinite(a) and np.isfinite(b) and a > 0.0 >= b:   # check if it crosses
        x0, x1 = u_sweep_val[i-1], u_sweep_val[i]
        U_flutter = x0 - a * (x1 - x0) / (b - a)             # linear interp of zero crossing
        break
U_flutter_ms = None if U_flutter is None else float(cal.ustar_to_ms(U_flutter))

# step 9: plot (dimensional: Hz vs m/s) with Fig. 6 overlay
fig, (ax_f, ax_d) = plt.subplots(2, 1, figsize=(6, 7))
ax_f.plot(u_ms, f_hz[:, 0], label="QS mode 1")
ax_f.plot(u_ms, f_hz[:, 1], label="QS mode 2")
ax_f.plot([0, 0], F_U0_HZ, "ks", ms=6, label="exp. U=0 (Fig. 6)")
ax_f.set_ylabel("frequency [Hz]")

ax_d.axhline(0.0, color="k", lw=0.8)
ax_d.plot(u_ms, damp[:, 0]); ax_d.plot(u_ms, damp[:, 1])
ax_d.set_ylabel(r"damping ratio $\zeta$"); ax_d.set_xlabel("airspeed $U$ [m/s]")

# experimental flutter / fold markers
for ax in (ax_f, ax_d):
    ax.axvline(U_FLUTTER_MS, color="r", ls="--", lw=1.0, label=f"exp. flutter {U_FLUTTER_MS} m/s")
    ax.axvline(U_FOLD_MS,    color="r", ls=":",  lw=1.0, label=f"exp. fold {U_FOLD_MS} m/s")

if U_flutter_ms is not None:
    for ax in (ax_f, ax_d):
        ax.axvline(U_flutter_ms, color="C3", ls="-.", lw=1.0)
    ax_f.set_title(f"QS flutter at U = {U_flutter_ms:.1f} m/s  (exp. {U_FLUTTER_MS})")
else:
    ax_f.set_title("QS: no flutter crossing — pitch mode unstable from low U")

ax_f.legend(fontsize=8)
plt.tight_layout()
plt.show()