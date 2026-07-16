"""Stage 2 (Peters finite-state inflow) -- calibrate beta to the LCO amplitude anchor.

SCOPE NARROWED, 16 July 2026. This script used to claim it located the LCO fold and
a subcritical window. It cannot, and both numbers were artifacts. Anything about the
branch at or below the Hopf point now comes from examples/stage2_LCO_continuation.py.
What remains here is the beta calibration and the stable branch above the Hopf.

--- what this file is still trusted for ---
    - calibrate beta ONCE to the near-flutter pitch amplitude (14 deg @ 13.5 m/s,
      Garcia Perez Fig. 9), then hold it fixed across every stage so that aero
      fidelity stays the only independent variable.
    - trace the stable LCO branch ABOVE the Hopf, where marching settles cleanly --
      verification controls A-D pass there (invariant to tolerance 1e-5..1e-10,
      to RK45/DOP853/Radau, and to initial kicks of 2..20 deg).

--- the trap, named ---
    The original docstring stated, correctly, that "near flutter the LCO settles so
    slowly that a fixed-window collapse test cannot separate slow growth from decay,
    so the fold comes out numerically unstable" -- and then reported "predicted fold
    ~13.0 m/s" as a result two paragraphs later. Both statements sat in this file at
    the same time. The diagnosis was right; the number was reported anyway.

    Mechanism, now fixed. _settle grew its averaging window until the two halves of
    the settled region agreed, then gave up at tau = 3500. The give-up branch was
    commented "flagged by caller" -- but the return signature carried no flag, so no
    caller could ever have flagged it. Worse, the loop also broke early whenever the
    amplitude fell under COLLAPSE_DEG, i.e. it declared convergence exactly where the
    window most needed to grow. _down_sweep then took the first speed whose
    (unconverged, still-decaying) amplitude fell under that threshold and called it
    the fold. That number was a stopwatch reading: the speed at which a transient
    happened to drop below a hard-coded threshold inside a hard-coded time budget.
    It was never a bifurcation.

    _settle now returns `converged` and callers honour it. The fold detection is
    GONE rather than flagged -- flagging is not enough, because the settling rate
    genuinely vanishes at a fold, so in finite time the quantity does not exist. No
    tau cap, tolerance or threshold repairs that. Shooting removes the question
    instead of tuning it: a limit cycle is a state that returns to itself,
    y(T; y0) = y0, solved directly by Newton with no settling test at all.

--- what the model actually does ---
    Per stage2_LCO_continuation.py: the Stage 2 branch is born AT the Hopf
    (13.1485 m/s, against the eigenvalue sweep's independent 13.148) and grows
    continuously above it. No fold, no bistable window. The model is SUPERCRITICAL;
    the rig is SUBcritical with a fold at 11.85 m/s. That is a disagreement about the
    KIND of bifurcation, not a gap in a number -- and a hardening cubic spring with
    attached-flow aero yields supercritical necessarily. It was never going to fold.

--- calibration discipline (unchanged; now the only claim this file makes) ---
    Amplitude is the anchor because it is well-conditioned: it settles cleanly above
    flutter and is monotone in beta, so brentq brackets. The fold was never an
    available anchor -- see above. Calibrated beta ~ 1.33 at N = 6.
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
U_HOPF_PRED  = 13.148         # this build's linear flutter (eigenvalue sweep)

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

    Returns (pitch_deg, pitch_sd_deg, plunge_mm, settled_struct_state, converged).
    The state is seeded into the next (lower) speed for continuation.

    converged=False means the window hit the tau cap without the two halves agreeing,
    so the amplitude returned is a TRANSIENT, not a settled orbit. Callers must not
    treat it as a branch point. This is the flag the old give-up branch claimed to
    raise and never carried -- see the module docstring.
    """
    converged = True
    while True:
        t_eval = np.linspace(0.0, tau_end, int(n_pts))
        sol = integrate(p, np.asarray(y0_struct, float), (0.0, tau_end),
                        aero=aero, U_star=U_star, t_eval=t_eval, rtol=RTOL)
        late = sol.t > 0.75 * tau_end
        t_late = sol.t[late]
        mid = t_late[0] + 0.5 * (t_late[-1] - t_late[0])
        a1, _ = _peak_amp(sol.y[1][late & (sol.t <= mid)])      # pitch, rad
        a2, sd = _peak_amp(sol.y[1][late & (sol.t > mid)])
        # NB no `a2 < COLLAPSE_DEG` short-circuit here: that declared convergence
        # exactly when the amplitude was small, i.e. near the fold, which is where
        # the window most needs to grow. A slow decay read as a settled collapse.
        if abs(a2 - a1) / max(a2, 1e-9) < tol:
            break
        tau_end *= 1.5
        n_pts *= 1.5
        if tau_end > 3500.0:
            converged = False                                   # give up, and SAY SO
            break
    plunge_mm = _peak_amp(sol.y[0][late & (sol.t > mid)])[0] * cal.b * 1e3   # xi = h/b
    return np.degrees(a2), np.degrees(sd), plunge_mm, sol.y[:4, -1].copy(), converged


def _amp_at(beta, U_ms):
    """Settled pitch amplitude (deg) at U_ms for the given beta (fresh aero).

    Hard-fails on an unconverged march: this feeds brentq, so an unconverged point
    would silently miscalibrate beta -- and beta is held fixed across every stage.
    """
    p = _section(beta)
    pitch, _, _, _, converged = _settle(p, PetersFinite(p, 0.0, N=N_INFLOW),
                                        cal.ms_to_ustar(U_ms), [0.0, KICK_RAD, 0.0, 0.0])
    if not converged:
        raise RuntimeError(
            f"beta anchor did not settle at {U_ms} m/s (beta={beta:.3f}). The anchor "
            f"must be a converged orbit; do not calibrate off a transient.")
    return pitch


def _down_sweep(beta, speeds_ms):
    """Trace the stable LCO branch downward from above the Hopf, re-seeding each
    step from the previous settled state. Returns rows[U, pitch, sd, plunge].

    NO fold detection. It used to return a fold bracket; that was the artifact --
    see the module docstring. A branch point that does not converge is dropped
    (NaN) and reported on stdout, never quietly kept. If the branch stops
    converging as you approach the Hopf, that is the method reaching its limit,
    not a bifurcation: use stage2_LCO_continuation.py.
    """
    p, aero = _section(beta), PetersFinite(_section(beta), 0.0, N=N_INFLOW)
    y0, rows = [0.0, KICK_RAD, 0.0, 0.0], []
    for U_ms in speeds_ms:
        pitch, sd, plunge, seed, converged = _settle(p, aero, cal.ms_to_ustar(U_ms), y0)
        if not converged:
            print(f"  [down] U = {U_ms:5.2f} m/s: UNRESOLVED at the tau cap -- dropped")
            rows.append([U_ms, np.nan, np.nan, np.nan])
        elif pitch < COLLAPSE_DEG:
            rows.append([U_ms, np.nan, np.nan, np.nan])         # decayed to zero
        else:
            rows.append([U_ms, pitch, sd, plunge])
            y0 = seed
    return np.array(rows)


def _up_sweep(beta, speeds_ms, tiny=np.radians(0.3)):
    """Cold-start each speed with a tiny IC (no re-seeding) to map the stable zero
    branch -- it stays ~0 below the Hopf. Unconverged points come back NaN.

    This no longer "exposes the bistable window": there is no window. The zero
    branch's stability below 13.148 is set by the eigenvalue sweep; this is only a
    cross-check of it.
    """
    p, aero = _section(beta), PetersFinite(_section(beta), 0.0, N=N_INFLOW)
    rows = []
    for U in speeds_ms:
        pitch, _, _, _, converged = _settle(p, aero, cal.ms_to_ustar(U),
                                            [0.0, tiny, 0.0, 0.0])
        if not converged:
            print(f"  [up]   U = {U:5.2f} m/s: UNRESOLVED at the tau cap -- dropped")
        rows.append([U, pitch if converged else np.nan])
    return np.array(rows)


# --- calibrate beta once, to the near-flutter amplitude (fold is then predicted) --
# amplitude decreases monotonically with beta, so brentq brackets cleanly
beta = brentq(lambda b: _amp_at(b, AMP_ANCHOR_MS) - AMP_TARGET_DEG, 0.2, 1.8, xtol=0.01)

# sweep from ABOVE the calibration anchor (13.5) downward, fine step. Points near
# the Hopf will stop converging -- that is expected and they are dropped, not read
# as a fold. The branch below the Hopf is stage2_LCO_continuation.py's job.
bif  = _down_sweep(beta, np.arange(14.0, 12.6, -0.1))
zero = _up_sweep(beta, np.arange(11.0, 13.4, 0.2))

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
    ax.axvline(U_HOPF_PRED,  color="0.5", ls="-.", lw=1.0, label=f"model Hopf {U_HOPF_PRED}")
    ax.set_xlabel("airspeed $U$ [m/s]")

ax_p.set_title(f"Stage 2 Peters LCO -- stable branch only (calibrated beta={beta:.2f})")
ax_p.legend(fontsize=8)
plt.tight_layout()
print(f"Ka_eff = {cal.Ka_eff:.2f} N.m/rad (nominal 11.53); "
      f"calibrated beta = {beta:.3f} at {AMP_TARGET_DEG:.0f} deg / {AMP_ANCHOR_MS} m/s")
live_n = int(np.sum(~np.isnan(bif[:, 1])))
print(f"Stage 2 Peters N={N_INFLOW}: {live_n}/{len(bif)} down-sweep points converged.")
print("No fold is reported: this method cannot measure one. See "
      "examples/stage2_LCO_continuation.py.")
plt.show()