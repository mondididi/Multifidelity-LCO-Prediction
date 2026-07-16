"""Stage 2 -- does the quadratic stiffness term flip the bifurcation subcritical?

stage2_LCO_continuation.py showed the Stage-2 model is SUPERcritical: the LCO branch
is born at the Hopf point (13.1485 m/s, matching the eigenvalue sweep's 13.148) and
grows continuously above it. The rig is SUBcritical -- fold at 11.85, bistable,
hysteretic. That is a disagreement about the KIND of bifurcation, not a gap in a
number, and a hardening cubic spring with attached-flow aero cannot produce a fold.

Lee et al. (1999) Prog. Aero. Sci. 35 Eq. 31 gives the general cubic spring as
M(a) = b0 + b1*a + b2*a^2 + b3*a^3. Stages 0-2 used the Woolston symmetric special
case -- b0 = b2 = 0 -- which is exact for a symmetric rig at zero trim. The Michigan
rig sits at 10 deg. This script asks what the missing b2 does.

Parameterisation. delta is the offset between the spring's neutral point and the
operating point; expanding a cubic spring about a point delta away regenerates
Lee's quadratic with b2 = 3*delta*b3. delta is NOT free -- it is bounded by the rig:
    delta = 0 deg   spring neutral AT the trim -> pure cubic -> the Stage-2 model
    delta = 10 deg  spring neutral at chord-zero, wing trimmed at 10 deg -> the
                    maximum physically available
The true delta is unknown (Garcia Perez state the 10 deg AoA and cite [26] for the
spring mechanism, but not where its neutral point sits -- that is Fig. 4). So sweep
the bracket and report the SENSITIVITY, not a fitted point.

Discipline. b1 = r_alpha_sq is held FIXED at the calibrated value, so the linear
system -- and therefore the flutter speed, 13.148 -- does not move with delta. delta
changes ONLY the nonlinear terms. Any change in bifurcation type is then attributable
to the nonlinearity alone. beta is RE-CALIBRATED at each delta to hold the 14 deg
amplitude anchor: one knob, one target, delta as the physical bracket. The fold
therefore remains a PREDICTION at every delta, not a fit.

What would count as fitting, and is not done here: choosing the delta that lands the
fold on 11.85 and reporting it as the answer. The output is the curve, not a point.
If the curve does cross 11.85 inside the bracket, the honest claim is a CONSISTENCY
statement -- "the observed fold is consistent with a neutral offset of X deg" -- not
a prediction.

Test. Sub vs super is decided at the Hopf point: a supercritical branch lives ABOVE
U_H, a subcritical one BELOW. So continuation only has to step a little past the
Hopf -- no need to trace each branch to its fold.

Run:  PYTHONPATH=src python examples/stage2_delta_sweep.py
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
from scipy.optimize import brentq

U_FLUTTER_MS = 13.19          # rig flutter (Garcia Perez et al. J063736)
U_FOLD_MS    = 11.85          # rig fold
U_FLUTTER_PRED = 13.148       # this build's linear flutter; delta-INVARIANT by design

N_INFLOW  = 6
U_ANCHOR  = 13.5              # amplitude anchor, above flutter where marching settles
AMP_TARGET_DEG = 14.0
DELTAS_DEG = [10.0, 6.67]   # bracket EDGE first: 10 deg is the decisive point.
                            # 0.0 and 3.33 already done -> both SUPERCRITICAL

cal = calibrate_michigan(zeta=structural_zeta())


def _section(beta, delta_rad):
    return TypicalSectionParameters(
        a=cal.a, x_alpha=cal.x_alpha, r_alpha_sq=cal.r_alpha_sq,
        omega_ratio=cal.omega_ratio, mu=cal.mu, beta=beta, delta=delta_rad,
        zeta_h=cal.zeta, zeta_alpha=cal.zeta)


def _settle_amp(beta, delta_rad, U_ms, tau_end=900.0):
    """Peak pitch [deg] late in a march. Only used to calibrate beta at U_ANCHOR,
    which is above flutter where controls A-D confirm marching settles cleanly."""
    p = _section(beta, delta_rad)
    aero = PetersFinite(p, 0.0, N=N_INFLOW)
    t_eval = np.linspace(0.0, tau_end, 40000)
    sol = integrate(p, [0.0, np.radians(10.0), 0.0, 0.0], (0.0, tau_end), aero=aero,
                    U_star=cal.ms_to_ustar(U_ms), t_eval=t_eval, rtol=1e-8)
    late = sol.y[1][sol.t > 0.85 * tau_end]
    pk, _ = find_peaks(late)
    return float(np.degrees(late[pk].mean())) if pk.size else 0.0


def _seed(beta, delta_rad, U_ms, tau_end=1200.0):
    p = _section(beta, delta_rad)
    aero = PetersFinite(p, 0.0, N=N_INFLOW)
    t_eval = np.linspace(0.0, tau_end, 60000)
    sol = integrate(p, [0.0, np.radians(10.0), 0.0, 0.0], (0.0, tau_end), aero=aero,
                    U_star=cal.ms_to_ustar(U_ms), t_eval=t_eval, rtol=1e-10)
    pk, _ = find_peaks(sol.y[1])
    if len(pk) < 2:
        raise RuntimeError("no orbit to seed from")
    return sol.y[:, pk[-1]].copy(), float(sol.t[pk[-1]] - sol.t[pk[-2]])


def _bail(tau, y, *a):
    """Terminal event: abandon a trial trajectory that has left the physical range.

    contique's dxmax caps the PREDICTOR step but not Newton's inner corrections, so
    an iterate can throw the state far from the orbit. At large pitch the cubic
    stiffness is violently stiff and RK45 grinds indefinitely -- the run appears to
    hang. Bailing out and returning a penalty keeps Newton moving instead.
    """
    return abs(y[1]) - np.radians(60.0)


_bail.terminal = True


def _shoot(x, lpf, p, aero, T_seed):
    y0, T = x[:10], x[10] * T_seed
    # y0 already outside the bail radius has no crossing for _bail to detect,
    # so the event never fires and the trial marches a full period at huge pitch.
    if (not (0.2 < x[10] < 5.0) or not np.all(np.isfinite(y0))
            or abs(y0[1]) > np.radians(60.0)):
        return np.full(11, 10.0)                    # push Newton back in-range
    sol = solve_ivp(structural_rhs, (0.0, T), y0, args=(p, aero, lpf),
                    method="RK45", rtol=1e-10, atol=1e-12, events=_bail)
    if sol.status == 1 or sol.t[-1] < 0.999 * T:    # bailed out
        return np.full(11, 10.0)
    return np.concatenate([sol.y[:, -1] - y0, [y0[3]]])


def _orbit(y0, T, p, aero, U_star):
    sol = solve_ivp(structural_rhs, (0.0, T), y0, args=(p, aero, U_star),
                    method="RK45", rtol=1e-10, atol=1e-12,
                    t_eval=np.linspace(0.0, T, 300))
    return np.degrees(sol.y[1].max()), np.degrees(sol.y[1].min())


# --- sweep ---------------------------------------------------------------------
print("=" * 76)
print("DELTA SWEEP -- does Lee's quadratic term flip the bifurcation subcritical?")
print("=" * 76)
print(f"  b1 = r_alpha_sq held FIXED -> flutter stays {U_FLUTTER_PRED} m/s for every delta")
print(f"  beta re-calibrated at each delta to {AMP_TARGET_DEG} deg @ {U_ANCHOR} m/s")
print(f"  rig: flutter {U_FLUTTER_MS}, fold {U_FOLD_MS}, subcritical\n")

out = {}
for d_deg in DELTAS_DEG:
    d = np.radians(d_deg)
    print(f"--- delta = {d_deg:5.2f} deg " + "-" * 45)
    try:
        beta = brentq(lambda b: _settle_amp(b, d, U_ANCHOR) - AMP_TARGET_DEG,
                      0.05, 6.0, xtol=0.01)
    except Exception as e:
        print(f"    calibration failed: {e}")
        continue
    b2_over_b3 = 3.0 * d
    print(f"    beta = {beta:.3f}   b2/b3 = 3*delta = {b2_over_b3:.4f} rad")

    y_s, T_s = _seed(beta, d, U_ANCHOR)
    p = _section(beta, d)
    aero = PetersFinite(p, 0.0, N=N_INFLOW)
    res = contique.solve(fun=_shoot, x0=np.append(y_s, 1.0),
                         lpf0=float(cal.ms_to_ustar(U_ANCHOR)),
                         args=(p, aero, T_s), dxmax=0.02, dlpfmax=0.005,
                         control0=(-1, -1), maxsteps=48, maxcycles=4,
                         maxiter=20, tol=1e-8, overshoot=1.05)
    rows = []
    for r in res:
        xs, lpf = r.x[:-1], float(r.x[-1])
        U = float(lpf * cal.b * cal.omega_alpha)
        amax, amin = _orbit(xs[:10], xs[10] * T_s, p, aero, lpf)
        rows.append([U, amax, amin])
    br = np.array(rows)
    out[d_deg] = (beta, br)

    # sub vs super: does the branch reach BELOW the Hopf point at finite amplitude?
    below = br[br[:, 0] < U_FLUTTER_PRED - 1e-4]
    finite = below[below[:, 1] > 1.0] if below.size else np.empty((0, 3))
    Umin, amin_at_turn = br[:, 0].min(), br[int(np.argmin(br[:, 0])), 1]
    print(f"    branch min U = {Umin:.4f} m/s at amplitude {amin_at_turn:.3f} deg")
    if finite.size:
        print(f"    => SUBCRITICAL: branch reaches {finite[:,0].min():.4f} m/s below "
              f"the Hopf at finite amplitude")
    else:
        print(f"    => SUPERCRITICAL: branch terminates at the Hopf, amplitude -> 0")

# --- report --------------------------------------------------------------------
print("\n" + "=" * 76)
print("SUMMARY")
print("=" * 76)
print(f"  {'delta':>7} | {'beta':>7} | {'b2/b3':>7} | {'min U':>8} | {'amp there':>9} | type")
print("  " + "-" * 68)
for d_deg, (beta, br) in out.items():
    Umin = br[:, 0].min()
    a_at = br[int(np.argmin(br[:, 0])), 1]
    typ = "SUBcritical" if (Umin < U_FLUTTER_PRED - 1e-4 and a_at > 1.0) else "SUPERcritical"
    print(f"  {d_deg:7.2f} | {beta:7.3f} | {3*np.radians(d_deg):7.4f} | {Umin:8.4f} | "
          f"{a_at:9.3f} | {typ}")
print(f"\n  rig fold {U_FOLD_MS} m/s; model flutter {U_FLUTTER_PRED} m/s (delta-invariant)")
np.save("/tmp/delta_sweep.npy", out, allow_pickle=True)

fig, ax = plt.subplots(figsize=(8.2, 5.0))
for d_deg, (beta, br) in out.items():
    ax.plot(br[:, 0], br[:, 1], "o-", ms=3, label=f"delta = {d_deg:.2f} deg (beta={beta:.2f})")
ax.axvline(U_FLUTTER_MS, color="r", ls="--", lw=1.2, label="rig flutter 13.19")
ax.axvline(U_FOLD_MS, color="r", ls=":", lw=1.4, label="rig fold 11.85")
ax.axvline(U_FLUTTER_PRED, color="0.5", ls="-.", lw=1.0, label="model flutter 13.148")
ax.set_xlabel("airspeed U [m/s]")
ax.set_ylabel("pitch LCO amplitude, max [deg]")
ax.set_title("Does Lee's quadratic term produce a fold?")
ax.legend(fontsize=8)
ax.grid(alpha=0.25)
plt.tight_layout()
plt.savefig("/tmp/delta_sweep.png", dpi=170)
print("  figure -> /tmp/delta_sweep.png")