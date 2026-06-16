# Michigan Section — Parameter Reference

Reference for `michigan_params.py`: where every number comes from, how the
derived ones are computed, the calibration derivation in full, and what each
test checks and why.

The Michigan rig is the García Pérez et al. (AIAA J) 2-DOF pitch–plunge wing —
our **experimental validation** target. Isogai Case A stays the numerical
**verification** benchmark; this section is calibrated to *measured* data.

---

## 1. The calibration principle (why we tune the structure at all)

At **U = 0 there is no aerodynamics**. The two measured frequencies (the U = 0
intercepts of Fig. 6, **5.3 / 6.2 Hz**) come purely from the mass and stiffness
matrices. So:

- They carry **zero information about aerodynamic fidelity** — every stage (QS,
  Peters, UVLM, CFD) gives the identical (zero) aero force at U = 0.
- They are a **structural fact** we *match*, not a result we *predict*.

We therefore calibrate the structure to them, which removes structural error as
a confound. The aero-fidelity result then lives entirely in the **U > 0**
comparison (frequency coalescence, flutter speed, LCO amplitude). We **never**
tune a structural parameter to match a U > 0 quantity — that would be tuning to
the result.

> Guitar analogy: tune to a reference pitch (U = 0 modes) *before* judging the
> playing (U > 0 aero). Nobody is graded on whether the guitar was in tune.

---

## 2. Parameter provenance (at a glance)

Three tiers: **rig facts** (read from the paper), **geometry-derived** (pure
algebra on the facts), **calibration-derived** (solved to hit the measured
modes). Plus two **input choices**.

| Symbol (code)            | Value      | Tier                 | Source / formula                |
|--------------------------|-----------:|----------------------|---------------------------------|
| `M_KG`                   | 2.4        | rig fact             | "total weight 2.4 kg"           |
| `S_ALPHA_KGM`            | 0.0285     | rig fact             | "static moment 0.0285 kg·m"     |
| `IA_BARE_KGM2`           | 3.401e-3   | rig fact*            | "pitch inertia 3.401e-3"        |
| `KH_NOM_NM`              | 1580       | rig fact*            | "1.58 N/mm" × 1000              |
| `KA_NOM_NMRAD`           | 11.53      | rig fact*            | "11.53 N·m/rad"                 |
| `B_M`                    | 0.1        | rig fact             | chord 200 mm ÷ 2                |
| `SPAN_M`                 | 0.44       | rig fact             | "wingspan 440 mm"               |
| `F_PLUNGE_HZ`            | 5.3        | measured (Fig. 6)    | U = 0 intercept, lower branch   |
| `F_PITCH_HZ`             | 6.2        | measured (Fig. 6)    | U = 0 intercept, upper branch   |
| `RHO_AIR`                | 1.225      | assumed              | std sea-level air (atmospheric) |
| `A_NONDIM`               | −0.5       | geometry-derived     | spar @ 25% chord (§4.1)         |
| `X_ALPHA`                | 0.11875    | geometry-derived     | `S_α/(m b)` (§4.2)              |
| `MU`                     | 141.7      | geometry-derived     | `m/(π ρ b² s)` (§4.3)          |
| `r_alpha_sq`             | 0.58267    | **calibration**      | solved → modes (§5)            |
| `omega_alpha`            | 35.80 rad/s| **calibration**      | solved → scale (§5)            |
| `Ia_eff / Kh_eff / Ka_eff`| —         | **calibration**      | back-out (§5.5)                 |
| `omega_ratio`            | 1.0        | input (pin)          | chooses family member (§5.4)    |
| `zeta`                   | 0.0        | input (TODO)         | recovery-rate back-out          |

**\* Kept but NOT used in the model.** `IA_BARE`, `KH_NOM`, `KA_NOM` are real rig
facts, but the structural model runs on the *calibrated effective* values
instead — bottom-up is wrong here (see §5.6 and §7). The nominals survive only
to print the ×-ratios as a sanity check.

---

## 3. Rig facts (the raw numbers)

From García Pérez et al. (NACA 0020, Bristol Prandtl tunnel):

```
m       = 2.4      kg          total wing mass
S_alpha = 0.0285   kg·m        static moment about the elastic axis
I_alpha = 3.401e-3 kg·m²       bottom-up pitch inertia about EA   (bare)
K_h     = 1580     N/m         plunge stiffness  (1.58 N/mm)      (nominal)
K_alpha = 11.53    N·m/rad     pitch stiffness                    (nominal)
chord   = 0.2      m   →  b = chord/2 = 0.1 m   (semichord)
span    = 0.44     m
spar    @ 25% chord  →  elastic axis location
```

Measured (read off Fig. 6 at U = 0):

```
f_plunge = 5.3 Hz   (lower / plunge-led branch)
f_pitch  = 6.2 Hz   (upper / pitch-led  branch)
```

---

## 4. Geometry-derived parameters

These are pure algebra on the rig facts. They are the parameters we **hold
fixed** during calibration (we trust geometry; we don't trust bottom-up
inertia/stiffness).

### 4.1 Elastic-axis location `a` (semichords from midchord, + toward TE)

The spar sits at 25% chord. We measure `a` from **midchord** (50% chord) in
units of semichord `b`:

```
a = (x_EA − x_mid) / b
  = (0.25 c − 0.50 c) / b
  = −0.25 c / b
  = −0.25 (2b) / b          (since c = 2b)
  = −0.5
```

Negative ⇒ EA is forward of midchord (toward the leading edge), as expected for
a quarter-chord spar.

### 4.2 Static imbalance `x_alpha`

The static moment is `S_α = m · d`, where `d` is the EA→CG distance. Nondimensionalising by `m b`:

```
x_alpha = d / b = S_alpha / (m b)
        = 0.0285 / (2.4 × 0.1)
        = 0.11875
```

Small ⇒ the CG is only ~0.12 semichords (≈6% chord) aft of the EA ⇒ weak
plunge–pitch coupling. (This matters in §5.4.)

### 4.3 Mass ratio `mu`

Wing mass vs. the mass of the air cylinder it sweeps (radius `b`, length `span`):

```
mu = m / (π ρ b² s)
   = 2.4 / (π × 1.225 × 0.1² × 0.44)
   = 141.7
```

Large `mu` ⇒ heavy wing relative to the air ⇒ aero is a weak perturbation on the
structure (consistent with low-speed flutter at modest U*).

---

## 5. The calibration (full derivation)

**Goal:** choose the *effective* structural parameters so the model's U = 0
coupled modes equal the measured 5.3 / 6.2 Hz, while holding the geometry
(`a`, `x_alpha`, `mu`) fixed.

### 5.1 The U = 0 structural eigenproblem (nondimensional)

State `q = [h/b, α]`, time nondimensionalised by `ω_α`. With

```
x  = x_alpha        wr = omega_ratio = ω_h / ω_α        r = r_alpha_sq
```

the mass and stiffness matrices are

```
M = [[ 1,   x  ],          K = [[ wr²,  0 ],
     [ x,   r  ]]               [ 0,    r ]]
```

The free-vibration eigenproblem `det(K − λ M) = 0`, with `λ = (ω/ω_α)²`:

```
det([[ wr² − λ,   −λ x      ],
     [ −λ x,      r − λ r    ]]) = 0

(wr² − λ)(r − λ r) − (λ x)² = 0
→  (r − x²) λ²  −  r (wr² + 1) λ  +  r wr²  =  0      ... (★)
```

So the two eigenvalues satisfy

```
λ₁ + λ₂ = r (wr² + 1) / (r − x²)        (sum)
λ₁ · λ₂ = r  wr²       / (r − x²)        (product)
```

### 5.2 Two facts from two frequencies

Each coupled frequency is `f_i = (ω_α/2π) √λ_i`, so `λ_i ∝ f_i²`. The two
measured frequencies give two independent pieces of information:

- a **ratio** (scale-free): `ρ = λ_lo/λ_hi = (f_lo/f_hi)²` → fixes `r_alpha_sq`
- an **absolute scale**: the actual size of either mode → fixes `omega_alpha`

This is the heart of it: **`r_alpha_sq` controls how close the two modes are;
`omega_alpha` controls where they sit on the Hz axis.**

### 5.3 Closed form for `r_alpha_sq`

Form the scale-free combination (depends only on the ratio):

```
T = (λ₁ + λ₂)² / (λ₁ λ₂) = ρ + 2 + 1/ρ
```

Substitute the sum/product from §5.1:

```
T = r (wr² + 1)² / [ (r − x²) wr² ]
```

Solve for `r`:

```
r_alpha_sq = T wr² x²  /  [ T wr² − (wr² + 1)² ]      ... 
```

Then solve the quadratic (★) numerically for `λ_lo, λ_hi`, and recover the scale:

```
ω_α = 2π f_lo / √λ_lo  =  2π f_hi / √λ_hi              ... 
```

(The two expressions agreeing is the calibration self-consistency check.)

### 5.4 The veering obstruction and the `omega_ratio` pin

There are **3** effective unknowns (`r`, `ω_α`, and the split between `ω_h`/`ω_α`)
but only **2** equations (the two modes). That leaves a **one-parameter family**.
We pin it with `omega_ratio = ω_h/ω_α`:

- `omega_ratio = 1.0` (default): `ω_h = ω_α` — the rig sits exactly at **veering**
  (uncoupled frequencies coincident). The 1.17 measured ratio is itself evidence
  for this.
- `omega_ratio ≈ 0.95`: nudges the lower mode to be genuinely **plunge-led**,
  matching the paper's branch labelling. Inertia rises slightly; scale barely
  moves.

**Why a clean 2-knob solve fails.** With `x_alpha` fixed at the geometric 0.119,
the two modes physically cannot get as close as 1.17 until the effective inertia
is ~4× bare — they're pinned near the closest-approach (veering) point. Any
attempt to hold `x_alpha` **and** a measured stiffness and solve for the other
two drives a parameter negative. The fix is to let **`ω_α` (the time scale)
float** as the third degree of freedom — which is exactly what this solver does.
The `denom ≤ 0` guard catches a pin so far from 1.0 that no
positive `r` exists.

### 5.5 Back-out of the effective dimensional structure

Pure definitions, inverted:

```
I_alpha,eff = r · m · b²
K_alpha,eff = I_alpha,eff · ω_α²          (from ω_α² = K_α / I_α)
K_h,eff     = m · (wr · ω_α)²             (from ω_h = wr ω_α, ω_h² = K_h/m)
```

These are record-only: they drive the ×-ratios and let the test suite rebuild
the dimensional eigenproblem independently of the nondim solve.

### 5.6 Result (omega_ratio = 1.0)

```
a           = −0.5         (geometry)
x_alpha     =  0.11875     (geometry)
r_alpha_sq  =  0.58267     (calibrated;  bottom-up was 0.14171)
mu          =  141.73      (geometry)
omega_ratio =  1.0         (veering pin)
omega_alpha =  35.80 rad/s (5.70 Hz)

effective vs bottom-up:  I_α 4.11×   K_h 1.95×   K_α 1.55×
```

The structure needs a big, broad correction — *all three* of inertia and both
stiffnesses shift up substantially. That is consistent with García Pérez
themselves using **data-driven identification** rather than bottom-up: the
as-built rig (spar, bearings, encoder, the parallel-spring cubic mechanism)
carries effective inertia/stiffness the component sum misses.

---

## 6. Dimensional scales & conversions

`omega_alpha` and `b` set the axes for overlaying experimental data:

```
f [Hz]   = (ω/ω_α) · ω_α/(2π)        # nondim freq → Hz      (omega_to_hz)
U [m/s]  = b · ω_α · U*              # nondim vel  → m/s      (ustar_to_ms)
U*       = U / (b · ω_α)             # m/s → nondim           (ms_to_ustar)

with b = 0.1 m, ω_α = 35.80 rad/s  →  U [m/s] = 3.58 · U*
```

Experimental landmarks mapped to U\*:

```
subcritical fold  11.85 m/s  →  U* ≈ 3.31
linear flutter    13.19 m/s  →  U* ≈ 3.69
```

So sweep U\* over ~0–4 to bracket the experiment.

---

## 7. Sensitivity — read Fig. 6 carefully

Because the rig sits at veering, `r_alpha_sq` is **hypersensitive** to the exact
U = 0 intercepts, while the velocity scale is **robust**:

| Fig-6 reading | required `I_eff` | U-scale |
|---------------|-----------------:|--------:|
| 5.3 / 6.2     | 4.1×             | 3.58    |
| 5.0 / 6.5     | 1.5×             | 3.52    |
| 4.8 / 6.8     | 0.9×             | 3.48    |
| 5.5 / 6.0     | 13×              | 3.60    |

**Action:** digitise the Fig. 6 intercepts before committing `r_alpha_sq`. The
flutter location in U\* (~3.7) holds regardless of the reading, so piece 3
(the VGBF sweep) can proceed in parallel.

---

## 8. Test suite — what each test checks and why

In `test_michigan_params.py`. Each test rebuilds the **dimensional** eigenproblem
from the effective structure, independently of the nondim solver, so a test
failure points at a real disagreement rather than re-using the solver's own math.

### `test_u0_modes_match_measured`
**What:** the calibrated section reproduces 5.3 / 6.2 Hz at U = 0 (abs 1e-3).
**Why:** the central contract. The calibration exists to hit these two targets;
if this fails, the solver or the back-out is wrong. Everything else is secondary
to this.

### `test_mode_ratio_matches`
**What:** the mode *ratio* `f_hi/f_lo` matches `6.2/5.3` (rel 1e-4).
**Why:** isolates the scale-free quantity — the one the calibration is most
sensitive to near veering. It separates failure modes: a wrong **placement**
(ω_α) leaves the ratio right; a wrong **separation** (r_alpha_sq) breaks it. Makes
a regression report diagnostic, not just red/green.

### `test_geometry_held_fixed`
**What:** `a = −0.5`, `x_alpha ≈ 0.11875`, `mu ≈ 141.7`.
**Why:** enforces the discipline boundary. Geometry must be an *input*, never
silently retuned by the calibration. If a refactor let the solver touch
`x_alpha` or `mu`, this catches it.

### `test_velocity_scale_robust`
**What:** for three plausible Fig-6 readings, the U-scale `b·ω_α` stays in
(3.3, 3.8) m/s per U\*.
**Why:** encodes the §7 finding as a guarantee — inertia is fragile but the
velocity scale (hence the flutter U\*) is robust. Locks the claim that the U > 0
comparison doesn't hinge on the exact intercept reading; flags it loudly if a
future change ever breaks that.

### `test_family_members_all_reproduce_modes` (parametrised wr = 0.95, 1.0, 1.05)
**What:** every in-family `omega_ratio` pin still hits 5.3 / 6.2 Hz.
**Why:** the pin chooses *which* member of the one-parameter family (how
plunge-/pitch-led the modes are), but a valid family means **every** member
reproduces the measured modes. Confirms the pin only redistributes mode
character, never breaks the mode match.

### `test_section_object_built` (`xfail`, non-strict)
**What:** once `michigan_section()` is wired to `mflco.TypicalSection`, it returns
a real section exposing `system_matrix`.
**Why:** marks the one remaining integration TODO. `xfail` (non-strict) keeps CI
green now and **flips to a genuine check** the moment the section is wired —
remove the marker then to make it load-bearing.

--