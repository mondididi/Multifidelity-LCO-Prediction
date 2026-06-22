"""Stage 2 (Peters) VGBF on Isogai Case A -- linear flutter-boundary verification.

Isogai Case A is the default TypicalSectionParameters (a=-2.0, x_alpha=1.8,
r_alpha_sq=3.48, mu=60). It is a CODE-verification target, not a quantitative one:
Case A's signature physics is the *transonic dip* (the flutter speed plunging near
M~0.8 due to shock dynamics), which PG-corrected linear theory cannot capture --
that gap is the Stage 4 CFD story. What Peters should get here is the U*=0 oracle
exactly and a sensible linear flutter boundary.

Two things to know vs the Michigan script:
  * needs a wide U* range -- with correct C(k) the wake lag stabilises Case A, so
    flutter sits high (U*~17 incompressible, ~13 at M=0.8). QS over-destabilises
    (no C(k)) and flutters spuriously low (U*~8-11); Peters removes that.
  * needs a little structural damping -- with zeta=0 the modes sit on the axis at
    U*=0 and QS/Peters push them unstable from U*~0+, so there is no + -> - cross.
  * the U*=0 oracle compares against the *damped* structural modes (zeta shifts
    omega_d = omega_n*sqrt(1-zeta^2)), so undamped_natural_frequencies is not the
    right reference when zeta>0.
"""
import numpy as np
import matplotlib.pyplot as plt

from mflco.model.params import TypicalSectionParameters
from mflco.aero.peters_finite import PetersFinite
from mflco.aero.quasi_steady import QuasiSteady
from mflco.model.analysis import modal_analysis, system_matrix

M_INF   = 0.0      # incompressible; try 0.8 for the (linear) compressible boundary
ZETA    = 0.01     # small structural damping so a flutter crossing exists
N_INFLOW = 6
U_MAX   = 25.0

p = TypicalSectionParameters(zeta_h=ZETA, zeta_alpha=ZETA)   # Isogai Case A
aero = PetersFinite(p, M_inf=M_INF, N=N_INFLOW)

u_sweep = np.linspace(0.0, U_MAX, 400)
freq = np.full((len(u_sweep), 2), np.nan)
damp = np.full((len(u_sweep), 2), np.nan)
for i, U in enumerate(u_sweep):
    try:
        freq[i], damp[i] = modal_analysis(p, 0.0, U, aero)
    except ValueError:
        continue


def flutter(u, mind):
    for i in range(1, len(u)):
        a, b = mind[i - 1], mind[i]
        if np.isfinite(a) and np.isfinite(b) and a > 0.0 >= b:
            return u[i - 1] - a * (u[i] - u[i - 1]) / (b - a)
    return None


U_flutter = flutter(u_sweep, np.nanmin(damp, axis=1))

# U*=0 oracle vs the DAMPED structural modes
Qs0 = system_matrix(p, 0.0, 0.0, None)
ev0 = np.linalg.eigvals(Qs0)
struct_wd = np.sort(np.abs(ev0[ev0.imag > 1e-9].imag))
f0, _ = modal_analysis(p, 0.0, 0.0,
                       PetersFinite(p, M_inf=M_INF, N=N_INFLOW, include_apparent_mass=False))
assert np.allclose(np.sort(f0), struct_wd, rtol=1e-4), (np.sort(f0), struct_wd)

fig, (ax_f, ax_d) = plt.subplots(2, 1, figsize=(6, 7))
ax_f.plot(u_sweep, freq[:, 0], label=f"Peters N={N_INFLOW} mode 1")
ax_f.plot(u_sweep, freq[:, 1], label=f"Peters N={N_INFLOW} mode 2")
ax_f.set_ylabel(r"frequency $\omega/\omega_\alpha$")
ax_d.axhline(0.0, color="k", lw=0.8)
ax_d.plot(u_sweep, damp[:, 0]); ax_d.plot(u_sweep, damp[:, 1])
ax_d.set_ylabel(r"damping ratio $\zeta$"); ax_d.set_xlabel(r"reduced velocity $U^*$")
if U_flutter is not None:
    for ax in (ax_f, ax_d):
        ax.axvline(U_flutter, color="C2", ls="-.", lw=1.0, label=f"flutter U*={U_flutter:.1f}")
    ax_f.set_title(f"Isogai Case A (M={M_INF}) -- Peters flutter at U* = {U_flutter:.2f}")
else:
    ax_f.set_title(f"Isogai Case A (M={M_INF}) -- no flutter in U* < {U_MAX}")
ax_f.legend(fontsize=8)
plt.tight_layout()

# QS comparison (the spurious low flutter from missing C(k))
def flutter_for(aero_obj):
    d = np.full((len(u_sweep), 2), np.nan)
    for i, U in enumerate(u_sweep):
        try: _, d[i] = modal_analysis(p, 0.0, U, aero_obj)
        except ValueError: pass
    return flutter(u_sweep, np.nanmin(d, axis=1))

print(f"U*=0 oracle (vs damped struct modes): {np.sort(f0)}  ==  {struct_wd}")
print(f"Isogai Case A (M={M_INF}, zeta={ZETA}):")
print(f"  QS     flutter U* = {flutter_for(QuasiSteady(p, M_INF))}   (spurious -- no C(k))")
print(f"  Peters flutter U* = {U_flutter}")

plt.show()