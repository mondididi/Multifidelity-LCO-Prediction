"""Stage 1 (quasi-steady) VGBF on Isogai Case A — verification.

Verifies the QS aero + eigenvalue machinery against the known benchmark:
sweep U*, build the system matrix, take eigenvalues, plot frequency & damping
vs U*. Section = Isogai Case A (defaults); aero = quasi-steady, M=0.
"""
# NOTE: with zeta=0 (canonical Isogai) the torsion mode is unstable from U*~=0
# no structural damping + no C(k) lag = QS over-destabilises
# so "no flutter crossing" here is correct, not a bug. 
# For a figure with a clean coalescence + crossing, rerun with TypicalSectionParameters(zeta_h=0.01, zeta_alpha=0.01).

from mflco.aero.quasi_steady import QuasiSteady
from mflco.model.params import TypicalSectionParameters
from mflco.model.analysis import modal_analysis, undamped_natural_frequencies
import numpy as np
import matplotlib.pyplot as plt

# initialise
p  = TypicalSectionParameters()              # beta = 0, everything dafault = isaogai, michigan requires using michigan properties, section
qs = QuasiSteady(p, M_inf=0.0)               # incompressible first; swap to 0.8 later, stage

u_sweep_val = np.linspace(0.0, 30.0, 60)     # non-dim U*, not m/s

freq = np.full((len(u_sweep_val), 2), np.nan)   # (n_speeds, 2 modes), shape as (60,2)
damp = np.full((len(u_sweep_val), 2), np.nan)

# steps 1-7: sweep, eigenvalues, reduce to (freq, zeta), store, skip coalescence
for i, U_star in enumerate(u_sweep_val):
    try:
        freq[i], damp[i] = modal_analysis(p, alpha_eq=0.0, U_star=U_star, aero=qs) #output as conjugate pairs
    except ValueError:        # past coalescence -> non-oscillatory; leave NaN
        continue

# step 10: U*=0 oracle (aero skipped at zero speed -> must equal structural freqs)
nat, _ = undamped_natural_frequencies(p)    #structural
assert np.allclose(np.sort(freq[0]), np.sort(nat)), (freq[0], nat)  #at U* = 0.0, it must be close to structural freqs

# step 8: flutter = lowest U* where any mode's damping crosses + -> -
min_damp = np.nanmin(damp, axis=1)
U_flutter = None
for i in range(1, len(u_sweep_val)):
    a, b = min_damp[i-1], min_damp[i]
    if np.isfinite(a) and np.isfinite(b) and a > 0.0 >= b:  #check if it crosses
        x0, x1 = u_sweep_val[i-1], u_sweep_val[i]
        U_flutter = x0 - a * (x1 - x0) / (b - a)     # linear interp of zero crossing
        break

# step 9: plot
fig, (ax_f, ax_d) = plt.subplots(2, 1, sharex=True, figsize=(6, 7))
ax_f.plot(u_sweep_val, freq[:, 0], label="mode 1")
ax_f.plot(u_sweep_val, freq[:, 1], label="mode 2")
ax_f.set_ylabel(r"$\omega\,/\,\omega_\alpha$"); ax_f.legend()

ax_d.axhline(0.0, color="k", lw=0.8)
ax_d.plot(u_sweep_val, damp[:, 0]); ax_d.plot(u_sweep_val, damp[:, 1])
ax_d.set_ylabel(r"damping ratio $\zeta$"); ax_d.set_xlabel(r"$U^*$")

if U_flutter is not None:
    for ax in (ax_f, ax_d):
        ax.axvline(U_flutter, color="r", ls="--", lw=0.8)
    ax_f.set_title(f"flutter at U* = {U_flutter:.3f}")
else:
    ax_f.set_title("no flutter crossing — least-stable mode never positive")

plt.tight_layout()
plt.show()