"""Stage 2 (Peters finite-state inflow) LCO bifurcation on the Michigan section.

Same calibrated structure as stage2_michigan_eigenvalue_sweep.py (U=0 modes pinned
to 5.3/6.2 Hz, flutter ~13.1 m/s), now driven nonlinear by the cubic pitch spring
K_alpha = r_alpha^2 (1 + beta*alpha^2) and time-marched to the settled limit cycle.
Forward integration traces the STABLE upper LCO branch by continuation and locates
the fold; the UNSTABLE branch between fold and flutter is unreachable by time
marching (Garcia Perez recover it with data-driven forecasting) -- so the scope here
is stable branch + fold + the zero/LCO bistability window.

beta is fixed across all stages (aero fidelity is the only variable), so it is
calibrated ONCE here, to the experimental near-flutter pitch amplitude (~14 deg,
Fig. 9). The fold is then a PREDICTION, not a fit.

Calibration note -- amplitude vs fold cannot both be anchored. Pinning amplitude is
fast and well-conditioned (settles cleanly above flutter; amplitude is monotone in
beta). Pinning the fold is not: near flutter the LCO settles so slowly that a
fixed-window collapse test cannot separate slow growth from decay, so the fold comes
out numerically unstable. Hence amplitude is the anchor, fold the prediction.

Result (calibrated beta ~ 1.33, N=6):
    near-flutter amplitude   14 deg @ 13.5 m/s   (calibrated to Fig. 9)
    predicted fold           ~13.0 m/s           (exp 11.85)
The amplitude scale matches, but the predicted subcritical window (~0.2 m/s) is far
narrower than the rig's (~1.34 m/s): attached-flow Peters captures the LCO mechanism
and amplitude but under-predicts the hysteresis width. That gap is the Stage-2
finding, separate from the 10 deg AoA on a NACA 0020 (a thick section near
separation that attached-flow aero cannot see).
"""
import numpy as np
import matplotlib.pyplot as plt

from mflco.model.params import TypicalSectionParameters
from mflco.model.michigan_params import calibrate_michigan, structural_zeta
from mflco.aero.peters_finite import PetersFinite
from mflco.model.solver import integrate
from scipy.signal import find_peaks
from scipy.optimize import brentq

# experimental landmarks (Garcia Perez et al. J063736, dimensional)
U_FLUTTER_MS = 13.19
U_FOLD_MS    = 11.85
AMP_ANCHOR_MS = 13.5          # settle-able point just above flutter for calibration
AMP_TARGET_DEG = 14.0         # near-flutter pitch amplitude plateau (Fig. 9)

N_INFLOW     = 6              # Peters inflow states (matches the eigenvalue sweep)
KICK_RAD     = np.radians(10.0)   # large IC to land on the LCO branch
COLLAPSE_DEG = 0.2           # pitch below this = decayed to the zero equilibrium
RTOL         = 1e-6          # sweep tolerance; re-verify a point at 1e-8 if needed

cal = calibrate_michigan(zeta=structural_zeta())


def _section(beta):
    """Michigan section with cubic pitch coefficient beta (section_from_params
    gives beta=0; the nonlinear runs need it set, as in test_peters_timemarch)."""
    return TypicalSectionParameters(
        a=cal.a, x_alpha=cal.x_alpha, r_alpha_sq=cal.r_alpha_sq,
        omega_ratio=cal.omega_ratio, mu=cal.mu, beta=beta,
        zeta_h=cal.zeta, zeta_alpha=cal.zeta)


def _peak_amp(sig):
    """Garcia Perez's amplitude: mean (and spread) of the positive local maxima
    of the settled window; (0, 0) if the response has decayed to zero."""
    pk, _ = find_peaks(sig)
    a = sig[pk]
    return (float(a.mean()), float(a.std())) if a.size else (0.0, 0.0)


def _settle(p, aero, U_star, y0_struct, tau_end=500.0, n_pts=10000, tol=0.02):
    """March to the settled LCO. y0_struct is the 4 structural ICs ONLY -- the
    solver appends the inflow states (passing a pre-augmented vector raises in
    aero_rhs). The settled window is halved and the two halves' peak-mean
    amplitudes compared; if they disagree by > tol the window is grown, since the
    LCO settles slowly near flutter and a fixed window reads a transient as the LCO.

    Returns (pitch_deg, pitch_sd_deg, plunge_mm, settled_struct_state); the last is
    seeded into the next (lower) speed for continuation.
    """
    while True:
        t_eval = np.linspace(0.0, tau_end, int(n_pts))
        sol = integrate(p, np.asarray(y0_struct, float), (0.0, tau_end),
                        aero=aero, U_star=U_star, t_eval=t_eval, rtol=RTOL)
        late = sol.t > 0.75 * tau_end
        t_late = sol.t[late]
        mid = t_late[0] + 0.5 * (t_late[-1] - t_late[0])
        a1, _ = _peak_amp(sol.y[1][late & (sol.t <= mid)])      # pitch, rad
        a2, sd = _peak_amp(sol.y[1][late & (sol.t > mid)])
        if a2 < np.radians(COLLAPSE_DEG) or abs(a2 - a1) / max(a2, 1e-9) < tol:
            break
        tau_end *= 1.5
        n_pts *= 1.5
        if tau_end > 3500.0:
            break                                               # give up; flagged by caller
    plunge_mm = _peak_amp(sol.y[0][late & (sol.t > mid)])[0] * cal.b * 1e3   # xi = h/b
    return np.degrees(a2), np.degrees(sd), plunge_mm, sol.y[:4, -1].copy()


def _amp_at(beta, U_ms):
    """Settled pitch amplitude (deg) at U_ms for the given beta (fresh aero)."""
    p = _section(beta)
    return _settle(p, PetersFinite(p, 0.0, N=N_INFLOW),
                   cal.ms_to_ustar(U_ms), [0.0, KICK_RAD, 0.0, 0.0])[0]


def _down_sweep(beta, speeds_ms):
    """Trace the stable upper branch from above flutter down through the fold,
    reusing one aero object and re-seeding each step from the previous settled
    state (continuation). Returns (rows[U, pitch, sd, plunge], fold_bracket)."""
    p, aero = _section(beta), PetersFinite(_section(beta), 0.0, N=N_INFLOW)
    y0, prev, fold, rows = [0.0, KICK_RAD, 0.0, 0.0], None, None, []
    for U_ms in speeds_ms:
        pitch, sd, plunge, seed = _settle(p, aero, cal.ms_to_ustar(U_ms), y0)
        if pitch < COLLAPSE_DEG:
            if fold is None and prev is not None:
                fold = (U_ms, prev)
            rows.append([U_ms, np.nan, np.nan, np.nan])
        else:
            rows.append([U_ms, pitch, sd, plunge])
            y0, prev = seed, U_ms
    return np.array(rows), fold


def _up_sweep(beta, speeds_ms, tiny=np.radians(0.3)):
    """Cold-start each speed with a tiny IC (no re-seeding) to map the stable zero
    branch -- stays ~0 until flutter, exposing the bistable window."""
    p, aero = _section(beta), PetersFinite(_section(beta), 0.0, N=N_INFLOW)
    return np.array([[U, _settle(p, aero, cal.ms_to_ustar(U), [0.0, tiny, 0.0, 0.0])[0]]
                     for U in speeds_ms])


# --- calibrate beta once, to the near-flutter amplitude (fold is then predicted) --
# amplitude decreases monotonically with beta, so brentq brackets cleanly
beta = brentq(lambda b: _amp_at(b, AMP_ANCHOR_MS) - AMP_TARGET_DEG, 0.2, 1.8, xtol=0.01)

# sweep from ABOVE the calibration anchor (13.5) down through the fold, fine step,
# so the whole branch is traced -- starting below 13.5 truncates most of the curve
bif, fold = _down_sweep(beta, np.arange(14.0, 12.6, -0.1))
zero = _up_sweep(beta, np.arange(11.0, 13.4, 0.2))
fold_ms = None if fold is None else 0.5 * (fold[0] + fold[1])

fig, (ax_p, ax_l) = plt.subplots(1, 2, figsize=(11, 4.5))
live = ~np.isnan(bif[:, 1])
ax_p.errorbar(bif[live, 0], bif[live, 1], yerr=bif[live, 2], fmt="o-", capsize=3,
              label=f"Peters N={N_INFLOW}, beta={beta:.2f}")
ax_p.plot(zero[:, 0], zero[:, 1], "-", color="0.6", label="stable zero branch")
ax_l.plot(bif[live, 0], bif[live, 3], "o-", label=f"Peters N={N_INFLOW}, beta={beta:.2f}")

ax_p.set_ylabel("pitch LCO amplitude [deg]")
ax_l.set_ylabel("plunge LCO amplitude [mm]")
for ax in (ax_l, ax_p):
    ax.axvline(U_FLUTTER_MS, color="r", ls="--", lw=1.0, label=f"exp. flutter {U_FLUTTER_MS}")
    ax.axvline(U_FOLD_MS,    color="r", ls=":",  lw=1.0, label=f"exp. fold {U_FOLD_MS}")
    if fold_ms is not None:
        ax.axvline(fold_ms, color="C2", ls="-.", lw=1.0, label=f"pred. fold {fold_ms:.2f}")
    ax.set_xlabel("airspeed $U$ [m/s]")

ax_p.set_title(f"Stage 2 Peters LCO -- stable branch only (calibrated beta={beta:.2f})")
ax_p.legend(fontsize=8)
plt.tight_layout()
print(f"Ka_eff = {cal.Ka_eff:.2f} N.m/rad (nominal 11.53); "
      f"calibrated beta = {beta:.3f} at {AMP_TARGET_DEG:.0f} deg / {AMP_ANCHOR_MS} m/s")
print(f"Stage 2 Peters N={N_INFLOW}: predicted fold "
      f"{'none in range' if fold_ms is None else f'{fold_ms:.2f} m/s'}  (exp {U_FOLD_MS})")
plt.show()