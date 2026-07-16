"""Stage 2 -- LCO branch by shooting + arclength continuation (contique).

Replaces the forward-marching sweep in stage2_michigan_LCO_amplitude.py. That sweep
cannot resolve the fold: near flutter the settling rate goes to zero, so a fixed
window cannot separate slow growth from slow decay, and _settle's give-up branch
returns unconverged transients as if they were branch points (see
stage2_LCO_verification.py, control E -- every point in the window came back
unresolved). No tau cap or tolerance fixes that; the quantity being measured
genuinely vanishes at the fold.

Shooting removes the question. A limit cycle is a state that returns to itself:

    y(T; y0) = y0

Solve that for (y0, T) by Newton -- no marching-until-settled, no settling
tolerance, no collapse threshold, no tau cap. One period of integration per
residual (tau ~ 21) instead of tau = 500..3500 per sweep point.

contique supplies the outer loop: ARCLENGTH continuation of g(x, lpf) = 0. This is
the part that matters. Stepping U and re-solving cannot pass a fold -- there
dU/ds = 0, the branch turns back on itself, and stepping in U walks off the end of
the curve. Arclength parameterises by distance ALONG the branch, which stays smooth
through the turn, so continuation rounds the fold and comes back up the UNSTABLE
branch -- the one this repo's docstring calls unreachable by time marching, the one
Garcia Perez recover with data-driven forecasting (Fig. 17, middle branch).

contique knows nothing about LCOs or ODEs: it follows g(x, lpf) = 0 and nothing
more. The physics is entirely in _shoot below. Division of labour:
    contique      -- arclength branch-following, fold turning
    _shoot        -- the periodic-orbit residual (the physics)
    _seed_march   -- one settled march to start from, at 13.5 m/s where marching
                     still works (controls A-D of the verification battery)

Unknowns x = [y0 (4 structural + 6 Peters inflow), T/T_seed] -- 11 total.
Equations = 10 periodicity residuals + 1 phase condition. The phase condition is
required: a periodic orbit is invariant under time-shift, so without pinning the
phase the Jacobian is singular. alpha_dot(0) = 0 pins the start at a pitch
extremum, which is also where Garcia Perez read amplitude.

The period is normalised by the seed period so every unknown is O(1) -- contique's
dxmax caps all unknowns with one number, and a raw T ~ 21 against states ~ 0.25
would make the arclength step size meaningless.

Note -- integrate() in solver.py hardcodes the inflow ICs to zero, so it cannot be
used here: shooting needs the FULL augmented state as unknowns, since the inflow
states must return to themselves too. structural_rhs takes the full state directly,
so solve_ivp is called on it here. No production change needed.

Recording max AND min per orbit, not a single amplitude (cf. recordPO in
dawbarton/AeroelasticCBC.jl): if a quadratic stiffness term is ever added the orbit
becomes asymmetric and one amplitude number is lossy.

Run:  PYTHONPATH=src python examples/stage2_LCO_continuation.py
"""
import numpy as np
import matplotlib.pyplot as plt
import contique

from mflco.model.params import TypicalSectionParameters
from mflco.model.michigan_params import calibrate_michigan, structural_zeta
from mflco.aero.peters_finite import PetersFinite
from mflco.model.eom import structural_rhs
from mflco.model.solver import integrate
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks

# experimental landmarks (Garcia Perez et al. J063736, dimensional)
U_FLUTTER_MS = 13.19
U_FOLD_MS    = 11.85
U_FLUTTER_PRED = 13.148       # this build's linear flutter (eigenvalue sweep)

N_INFLOW  = 6                 # matches the eigenvalue sweep and the LCO sweep
BETA_CAL  = 1.326             # calibrated ONCE to 14 deg @ 13.5 m/s in
                              # stage2_michigan_LCO_amplitude.py; held fixed here.
                              # Re-run that script's brentq if the structure changes.
U_SEED_MS = 13.5              # above flutter, where marching still settles cleanly
SHOOT_RTOL = 1e-10            # tight: the residual IS the accuracy of the orbit

cal = calibrate_michigan(zeta=structural_zeta())


def _section(beta):
    """Michigan section with cubic pitch coefficient beta (section_from_params
    gives beta=0; the nonlinear runs need it set, as in test_peters_timemarch)."""
    return TypicalSectionParameters(
        a=cal.a, x_alpha=cal.x_alpha, r_alpha_sq=cal.r_alpha_sq,
        omega_ratio=cal.omega_ratio, mu=cal.mu, beta=beta,
        zeta_h=cal.zeta, zeta_alpha=cal.zeta)


def _seed_march(beta, U_ms, tau_end=1200.0):
    """One settled march to seed the continuation. Returns (y0_full, T).

    Uses integrate() -- legitimate here because this is a cold start from rest, so
    zero inflow ICs are correct, and 13.5 m/s is above flutter where controls A-D
    of the verification battery confirm the march settles cleanly. The seed only
    has to be inside Newton's basin; the shooting solve makes it exact.

    y0 is taken at the LAST pitch peak (where alpha_dot ~ 0, matching the phase
    condition) and T from the spacing of the final two peaks.
    """
    p = _section(beta)
    aero = PetersFinite(p, 0.0, N=N_INFLOW)
    t_eval = np.linspace(0.0, tau_end, 60000)
    sol = integrate(p, [0.0, np.radians(10.0), 0.0, 0.0], (0.0, tau_end), aero=aero,
                    U_star=cal.ms_to_ustar(U_ms), t_eval=t_eval, rtol=1e-10)
    pk, _ = find_peaks(sol.y[1])
    if len(pk) < 2:
        raise RuntimeError("no LCO to seed from -- march did not reach an orbit")
    T = float(sol.t[pk[-1]] - sol.t[pk[-2]])
    return sol.y[:, pk[-1]].copy(), T


def _shoot(x, lpf, p, aero, T_seed):
    """Periodic-orbit residual. contique drives x and lpf; lpf IS U_star.

    Returns 11 equations for 11 unknowns:
        y(T; y0) - y0   (10)  the orbit closes on itself
        alpha_dot(0)    (1)   phase condition -- start at a pitch extremum
    """
    y0 = x[:10]
    T = x[10] * T_seed
    sol = solve_ivp(structural_rhs, (0.0, T), y0, args=(p, aero, lpf),
                    method="RK45", rtol=SHOOT_RTOL, atol=1e-12)
    return np.concatenate([sol.y[:, -1] - y0, [y0[3]]])


def _orbit_minmax(y0, T, p, aero, U_star):
    """Re-integrate one period and record max/min pitch [deg] and period.

    max and min separately, not a single amplitude (cf. recordPO in
    AeroelasticCBC.jl) -- an asymmetric orbit is not described by one number.
    """
    sol = solve_ivp(structural_rhs, (0.0, T), y0, args=(p, aero, U_star),
                    method="RK45", rtol=1e-10, atol=1e-12,
                    t_eval=np.linspace(0.0, T, 400))
    return np.degrees(sol.y[1].max()), np.degrees(sol.y[1].min())


# --- seed the branch from one settled march at 13.5 m/s -----------------------
print(f"Ka_eff = {cal.Ka_eff:.2f} N.m/rad; beta = {BETA_CAL} (fixed); N = {N_INFLOW}")
print(f"Seeding from a settled march at {U_SEED_MS} m/s ...")
p = _section(BETA_CAL)
aero = PetersFinite(p, 0.0, N=N_INFLOW)
y_seed, T_seed = _seed_march(BETA_CAL, U_SEED_MS)
U_star_0 = float(cal.ms_to_ustar(U_SEED_MS))
print(f"  seed: pitch {np.degrees(y_seed[1]):.3f} deg, period tau = {T_seed:.3f}, "
      f"U* = {U_star_0:.4f}")

r0 = _shoot(np.append(y_seed, 1.0), U_star_0, p, aero, T_seed)
print(f"  seed residual |r| = {np.linalg.norm(r0):.3e}  (Newton will drive this to 0)")

# --- continue DOWN in U*, around the fold, up the unstable branch -------------
# control0 = (-1, -1): control the LAST extended unknown (lpf) in the NEGATIVE
# direction. contique switches the control component automatically at the fold --
# where lpf stops being a valid parameter, it hands control to whichever unknown is
# moving fastest. That switch IS how the branch turns.
print("Continuing (arclength) ...")
res = contique.solve(
    fun=_shoot,
    x0=np.append(y_seed, 1.0),
    lpf0=U_star_0,
    args=(p, aero, T_seed),
    dxmax=0.02,
    dlpfmax=0.004,            # ~0.015 m/s per step near the fold
    control0=(-1, -1),
    maxsteps=40,
    maxcycles=4,
    maxiter=20,
    tol=1e-8,
    overshoot=1.05,
)

# --- unpack: r.x is the extended vector [x (11), lpf] ------------------------
rows = []
for r in res:
    xs, lpf = r.x[:-1], float(r.x[-1])
    U_ms = float(lpf * cal.b * cal.omega_alpha)
    a_max, a_min = _orbit_minmax(xs[:10], xs[10] * T_seed, p, aero, lpf)
    rows.append([U_ms, a_max, a_min, xs[10] * T_seed])
br = np.array(rows)
np.savetxt("/tmp/branch.txt", br, header="U_ms  alpha_max_deg  alpha_min_deg  period_tau")
print("\n   U [m/s] | max [deg] | min [deg] | period")
for row in br:
    print(f"  {row[0]:8.4f} | {row[1]:9.4f} | {row[2]:9.4f} | {row[3]:7.3f}")

print(f"  {len(br)} branch points, U from {br[:, 0].min():.3f} to {br[:, 0].max():.3f} m/s")
turn = int(np.argmin(br[:, 0]))
print(f"  minimum U on the branch = {br[turn, 0]:.3f} m/s at step {turn}"
      f"  (experimental fold {U_FOLD_MS})")
if turn not in (0, len(br) - 1):
    print("  -> the branch TURNED: arclength rounded the fold and the points after"
          " the turn are the UNSTABLE branch.")
else:
    print("  -> no turn found in range; the branch ran to the end of the sweep.")

fig, (ax_b, ax_t) = plt.subplots(1, 2, figsize=(11, 4.5))
ax_b.plot(br[:turn + 1, 0], br[:turn + 1, 1], "o-", ms=3, label="stable branch (max)")
ax_b.plot(br[:turn + 1, 0], br[:turn + 1, 2], "o-", ms=3, color="C0", alpha=0.4,
          label="stable branch (min)")
if turn not in (0, len(br) - 1):
    ax_b.plot(br[turn:, 0], br[turn:, 1], "s--", ms=3, color="C3",
              label="unstable branch (max)")
ax_t.plot(br[:, 0], br[:, 3], "o-", ms=3, label="period")

ax_b.set_ylabel("pitch [deg]")
ax_t.set_ylabel(r"period $\tau$")
for ax in (ax_b, ax_t):
    ax.axvline(U_FLUTTER_MS, color="r", ls="--", lw=1.0, label=f"exp. flutter {U_FLUTTER_MS}")
    ax.axvline(U_FOLD_MS, color="r", ls=":", lw=1.0, label=f"exp. fold {U_FOLD_MS}")
    ax.axvline(U_FLUTTER_PRED, color="0.5", ls="-.", lw=1.0, label="pred. flutter")
    ax.set_xlabel("airspeed $U$ [m/s]")
ax_b.set_title(f"Stage 2 Peters LCO -- shooting + arclength (beta={BETA_CAL})")
ax_b.legend(fontsize=8)
plt.tight_layout()
plt.savefig("/tmp/branch.png", dpi=170)   # savefig, not show: this runs backgrounded