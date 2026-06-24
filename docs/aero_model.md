# Aerodynamic models

This document specifies the conventions, theory, and implementation contracts for the aerodynamic models in stages 1–4 of the multifidelity LCO project. Section 1 (conventions) is shared by all stages. Section 2 introduces the general unsteady aerodynamic problem that every model is approximating. Section 3 covers Stage 1 (quasi-steady with Prandtl-Glauert). Stages 2–4 will be added as they are implemented.

---

## 1. Conventions and notation

### 1.1 Geometry

The typical section is a 2D rigid airfoil of chord 2b, where b is the **semi-chord**. All geometric distances are measured in semi-chords from the **mid-chord**, with positive measured toward the trailing edge.

- **a** (also written **a_h** following Lee): non-dimensional position of the **elastic axis (EA)** from mid-chord. Negative if forward of mid-chord. For Isogai Case A, a = −2 (EA forward of the leading edge — an idealised configuration).
- **x_α**: non-dimensional position of the **centre of gravity (CG)** behind the EA. Positive if CG is aft of EA. For Isogai Case A, x_α = 1.8.
- Quarter-chord is at −0.5 (the aerodynamic centre for thin-airfoil theory).
- 3/4-chord (the "control point" for unsteady aero) is at +0.5; distance from EA is (½ − a) · b.

### 1.2 Degrees of freedom

- **h**: plunge displacement of the EA, **positive downward** (Lee 1999 convention).
- **α**: pitch angle, **positive nose-up**.
- **ξ ≡ h/b**: non-dimensional plunge.

State vector: q = [ξ, α]ᵀ, with q̇ = [ξ̇, α̇]ᵀ and q̈ = [ξ̈, α̈]ᵀ.

### 1.3 Time convention

This code uses

> τ = ω_α · t

where ω_α is the uncoupled pitch natural frequency. Primes denote derivatives with respect to this τ:

> ξ' = dξ/dτ,  ξ'' = d²ξ/dτ²

**Lee (1999) uses a different convention**: τ_Lee = U·t / b. These are related by

> τ_Lee = U* · τ

where U* = U / (b · ω_α) is the non-dimensional freestream speed. Consequently:

> ξ'_Lee = ξ' / U*,  ξ''_Lee = ξ'' / U*²

with the analogous relations for α. **Every expression sourced from Lee must be converted before use in code.** The unconverted Lee form is preserved in the derivations below for traceability; the final code-ready expressions are in this codebase's convention. See §1.7 for systematic procedure.

**Note: Theodorsen uses Frequency-domain via reduced frequency k**

### 1.4 Non-dimensional parameters

- **U*** = U / (b · ω_α): non-dimensional freestream speed. (U* - U star)
- **μ** = m / (π · ρ · b²): mass ratio (airfoil mass per unit span / mass of the cylinder of air with diameter 2b).
- **r_α²** = I_α / (m · b²): non-dimensional radius of gyration about the EA.
- **ω̄** = ω_h / ω_α: ratio of uncoupled plunge to pitch natural frequencies.
- **M_∞**: freestream Mach number.

### 1.5 Aerodynamic force outputs

Each aero model returns two scalar quantities, Q_ξ and Q_α, which enter the equations of motion on the right-hand side:

> M · q̈ + C · q̇ + K(α) · q = Q_aero,  Q_aero = [Q_ξ, Q_α]ᵀ

These are obtained from the physical lift coefficient C_L (positive upward) and moment coefficient about the EA C_M (positive nose-up) via

> Q_ξ = −(U*² / πμ) · C_L
> Q_α = +(2 · U*² / πμ) · C_M

The minus sign in Q_ξ reflects the sign-convention mismatch: C_L is positive upward, while ξ is positive downward, so positive lift produces a *negative* generalized force on the ξ coordinate.

### 1.6 Derivation of the generalized force prefactors

The formulas Q_ξ = −(U*²/πμ)·C_L and Q_α = +(2·U*²/πμ)·C_M used in §1.5 are not arbitrary. They follow from nondimensionalising the dimensional Newton's equations of motion for the typical section. This subsection traces the algebra. It applies to every aero model in this codebase — the prefactors are *structural*, not aerodynamic, and only the recipe for C_L and C_M changes from stage to stage.

#### Starting from the dimensional EOM

For the 2-DOF airfoil section in dimensional form (Lee 1999 Eqs. 8–9), Newton's law and angular-momentum balance give:

> m · ḧ + S_α · α̈ + C_h · ḣ + K_h · h = −L
> S_α · ḧ + I_α · α̈ + C_α · α̇ + K_α · α = M_EA

where m is mass per unit span, S_α = m·x_α·b is the static moment about the EA, I_α is the moment of inertia about the EA, C_h and C_α are damping constants, K_h and K_α are spring constants, L is lift (positive upward), and M_EA is moment about the EA (positive nose-up).

The minus sign on L appears because lift is positive *upward* while h is positive *downward*; positive lift produces negative h-acceleration.

#### Nondimensionalising the plunge equation

Divide the plunge equation through by m · b · ω_α². On the left, this produces ξ'' (since h/b = ξ and double-differentiation in τ = ω_α t introduces a 1/ω_α² factor), the static-moment coupling term, and structural damping/stiffness terms — all consistent with Stage 0's nondimensional EOM.

On the right, the aerodynamic term becomes:

> −L / (m · b · ω_α²)

Substitute the standard lift-coefficient definition L = ρ · U² · b · C_L (using full chord c = 2b and absorbing ½ into convention as Lee does), and substitute ρ = m / (π · μ · b²) from the definition μ = m / (π·ρ·b²):

> −L / (m · b · ω_α²)
> = −[m/(π·μ·b²)] · U² · b · C_L / (m · b · ω_α²)
> = −U² / (π · μ · b² · ω_α²) · C_L
> = −[U / (b·ω_α)]² · (1/(π·μ)) · C_L
> = **−(U*² / π·μ) · C_L**

The grouping U/(b·ω_α) is exactly U* (§1.4), so the prefactor collapses to U*²/(π·μ).

#### Nondimensionalising the pitch equation

Divide the pitch equation by I_α · ω_α². The right-hand side becomes:

> M_EA / (I_α · ω_α²)

Substitute the moment-coefficient definition M_EA = ρ · U² · b² · C_M · 2 (the factor of 2 reflects Lee's convention; some references absorb it into C_M itself — check the source). Substitute ρ = m/(π·μ·b²) and I_α = m · r_α² · b²:

> M_EA / (I_α · ω_α²)
> = 2 · [m/(π·μ·b²)] · U² · b² · C_M / (m · r_α² · b² · ω_α²)
> = 2 · U² / (π · μ · r_α² · b² · ω_α²) · C_M
> = **+(2·U*² / (π·μ·r_α²)) · C_M**

Depending on whether the r_α² is kept in the prefactor or moved into the EOM's left-hand side (as Lee does in his Eq. 11), the final form may absorb r_α² differently. In this codebase the convention used in §1.5 places r_α² in the structural mass matrix M[1,1] = r_α², so the Q_α prefactor is written as +(2·U*² / πμ)·C_M and the r_α² scaling is implicit through the EOM's mass matrix.

#### Summary

The two prefactors come from three pieces of dimensional analysis:

1. **U*²** — from rewriting U² in terms of (U/(b·ω_α))² (the nondimensional dynamic-pressure scaling).
2. **1/π** — from the definition of the mass ratio μ, which has a π in its denominator by convention.
3. **1/μ** — from the ratio of (fluid mass on scale b) to (airfoil mass m).

Each stage in this project will substitute its own expression for C_L and C_M into these prefactors. The prefactors themselves are model-independent.

#### Verification check

If C_L → 2π·α (steady thin-airfoil) and the airfoil is held at constant α with no motion, then Q_ξ → −(U*²/π·μ)·2π·α = −2·U*²·α/μ. This is the steady aerodynamic load on the plunge DOF. Equating it to a hypothetical equilibrium plunge spring force gives the divergence condition — sanity-check at the end of Stage 1 implementation.

### 1.7 Sign conventions, summarised

| Quantity      | Positive direction      |
|---------------|-------------------------|
| ξ = h/b       | downward                |
| α             | nose-up                 |
| C_L           | upward                  |
| C_M (about EA)| nose-up                 |
| Q_ξ           | downward (same as ξ)    |
| Q_α           | nose-up (same as α)     |

### 1.8 Working with other papers' conventions

Aeroelasticity literature is split between two time conventions:

- **Aerodynamics convention**: τ = U·t / b (chord-traversal time). Used by Theodorsen 1935, Fung 1955, Bisplinghoff 1955, Lee 1999, Jones 1940, Katz & Plotkin (UVLM). Natural for wake-development equations; reduced frequency k = ωb/U appears directly.
- **Structural convention**: τ = ω_n · t (vibration cycles), with ω_n a reference structural frequency. Used in classical structural dynamics texts.

This codebase uses the structural convention with ω_n = ω_α. Rationale:

1. **Structural baseline is independent of flow speed.** At U* = 0 (wind-off), the EOM reduces to a clean structural eigenvalue problem with no U* appearing — matching the Stage 0 verification plot.
2. **Stability analysis is transparent.** Eigenvalues of the assembled system can be watched as U* varies; U* = 0 recovers structural eigenvalues, and aerodynamic modification is identifiable term-by-term.
3. **Continuity with existing tests and `undamped_natural_frequencies`.** The structural model and its 11 verification tests were built in this convention; switching would require rewriting and re-verifying without conceptual benefit.
If the model is damped, then use 
**`system_matrix`** equa.

#### Procedure for porting expressions from a source paper

When importing equations from any reference into this codebase:

1. **Identify the source convention.** Look for explicit definitions of τ, U*, the velocity scale, and the length scale. Often stated near the EOM, sometimes only in a nomenclature table.
2. **Build a conversion table** between source primes and code primes (see Lee→Mond table below for the template).
3. **Substitute mechanically** into the source expression. This is straightforward but error-prone; do it on paper, then verify.
4. **Sanity-check at a known limit.** Set all rates to zero and confirm the steady result. Or set U* → 0 and confirm the result vanishes (for purely aerodynamic terms). Or compare against a published numerical result if one exists.

#### Lee (1999) → this codebase

| Source (Lee, τ_Lee = Ut/b) | This codebase (τ = ω_α · t) |
|---------------------------|------------------------------|
| τ                         | U* · τ                       |
| ξ'_Lee                    | ξ' / U*                      |
| ξ''_Lee                   | ξ'' / U*²                    |
| α'_Lee                    | α' / U*                      |
| α''_Lee                   | α'' / U*²                    |
| k = ωb/U                  | unchanged                    |
| C(k), φ(τ_Lee)            | C(k) unchanged; φ receives U*·τ as argument |
| U*, μ, r_α², ω̄, a, x_α    | unchanged                    |

**Note**: k is a property of the motion (ratio of oscillation rate to chord-traversal rate), not a convention choice — it has the same numerical value regardless of which τ is used.

#### Per-stage convention statement

Each stage's section in this document (and the corresponding class's docstring) should open with a statement of the form:

> "Expressions in this section are derived from [Source] Eqs. [X–Y] in [source convention]. The conversions in §1.7 have been applied. The implementation uses the code convention τ = ω_α · t throughout."

This makes the provenance of every expression traceable and prevents silent convention drift across stages.

---

## 2. The general unsteady aerodynamic problem

All four aero models in this project are approximations of the same underlying physics: the lift and moment on a thin airfoil undergoing pitch and plunge motion in a uniform freestream. Understanding what each model is an approximation *of* clarifies what it gets right and what it misses.

### 2.1 The four contributions to lift

For an airfoil pitching at rate α̇ and plunging at rate ḣ in a freestream of speed U, the lift has four physical contributions:

1. **Geometric angle of attack (α)** — the static lift from being tilted relative to the wind.
2. **Plunge-induced effective angle (ḣ / U)** — vertical motion changes the apparent wind direction by an angle ḣ/U in the airfoil frame.
3. **Pitch-rate effective angle ((½ − a) · b · α̇ / U)** — rotation makes different points along the chord move at different vertical velocities. Evaluating at the 3/4-chord (the kinematic control point) gives this contribution.
4. **Apparent mass** — accelerating the airfoil also accelerates the surrounding air, producing lift proportional to ḧ and α̈. This is a non-circulatory contribution (no bound vortex involved).

The first three combine into the **effective angle of attack**:

> α_eff = α + ḣ/U + (½ − a) · b · α̇ / U

In non-dimensional Lee-convention form (since ḣ/U = ξ'_Lee and b·α̇/U = α'_Lee):

> α_eff = α + ξ'_Lee + (½ − a) · α'_Lee

### 2.2 Theodorsen's function and the wake

Steady thin-airfoil theory says C_L = 2π · α. For an oscillating airfoil, this is modified by **Theodorsen's function** C(k):

> C_L_circ = 2π · C(k) · α_eff

where k = ω · b / U is the **reduced frequency**. Physically, when the airfoil's circulation changes, equal and opposite vorticity is shed into the wake by Kelvin's theorem. The shed vorticity drifts downstream but continues to induce velocities back on the airfoil. C(k) encodes the resulting:

- **amplitude attenuation**: |C(k)| < 1
- **phase lag**: arg C(k) < 0

Limits:

- C(0) = 1 (steady flow — wake fully dissipated).
- C(∞) = ½ (infinitely fast oscillation — wake right next to airfoil).

The phase lag in C(k) is what makes classical flutter possible. Without it, motion and aerodynamic force are in phase, and the energy exchange that drives flutter doesn't occur.

### 2.3 The fidelity ladder

Each stage in this project makes a different approximation of C(k):

| Stage | Aero model                            | Treatment of C(k)                              |
|-------|-------------------------------------  |------------------------------------------------|
| 1     | Quasi-steady + Prandtl-Glauert        | C(k) ≡ 1                                       |
| 2     | Peters finite-state inflow            | N induced-flow states; first-order wake lag (C(k) phase) |
| 3     | UVLM                                  | Discrete vortex tracking, exact within potential-flow assumptions |
| 4     | Unsteady Euler / intermediate CFD     | Full compressible inviscid flow                |

Stages 1–3 are all *inviscid, incompressible* in their raw form; Stage 1 adds compressibility via Prandtl-Glauert; Stage 4 handles compressibility natively.

---

## 3. Stage 1: quasi-steady aerodynamics with Prandtl-Glauert correction

### 3.1 Assumption

Quasi-steady aerodynamics is the limit φ(τ) → 1 of Lee's Eq. 15 — equivalently, C(k) ≡ 1. The wake is assumed to respond *instantaneously* to changes in airfoil motion. There is no memory, no phase lag, no aero state.

This is the cheapest possible aero model that has any unsteadiness at all. It still uses ḣ and α̇ to construct the effective angle of attack (a *kinematic* effect, not a wake effect), but it ignores the dynamic memory effect of the wake.

### 3.2 Derivation from Lee Eq. 15–16

Lee (1999), Eqs. 15–16, gives the full unsteady lift and moment in Lee's time convention as:

> C_L_Lee = π · (ξ''_L − a · α''_L + α'_L)
>          + 2π · {α(0) + ξ'_L(0) + (½ − a) · α'_L(0)} · φ(τ_L)
>          + 2π · ∫₀^τ_L  φ(τ_L − σ) · (dα_eff/dσ) dσ

> C_M_Lee = π · (½ + a) · {α(0) + ξ'_L(0) + (½ − a) · α'_L(0)} · φ(τ_L)
>          + π · (½ + a) · ∫₀^τ_L φ(τ_L − σ) · (dα_eff/dσ) dσ
>          + (π/2) · a · (ξ''_L − a · α''_L)
>          − (π/2) · (½ − a) · α'_L
>          − (π/16) · α''_L

With φ(τ) → 1, the indicial-IC term and the Duhamel integral combine, by the fundamental theorem of calculus, into simply 2π · α_eff(τ_L). The full QS lift coefficient in Lee convention is:

> **C_L_QS_Lee = π · (ξ''_L + α'_L − a · α''_L) + 2π · [α + ξ'_L + (½ − a) · α'_L]**

Similarly for the moment:

> **C_M_QS_Lee = π · (½ + a) · [α + ξ'_L + (½ − a) · α'_L]**
>              **+ (π/2) · a · (ξ''_L − a · α''_L)**
>              **− (π/2) · (½ − a) · α'_L**
>              **− (π/16) · α''_L**

The first term of each expression is **circulatory** lift/moment with C(k) = 1. The remainder is **non-circulatory** (apparent mass) contribution.

### 3.3 Conversion to code convention

Substituting ξ'_L = ξ' / U* and ξ''_L = ξ'' / U*² (and analogously for α), and then forming **Q_ξ = −U*² · C_L / (πμ)**  **Q_α = +(2·U*² / πμ) · C_M**
 (and analogously for Q_α), the **generalized forces** in this codebase's time convention become:  

> **Q_ξ = − (1/μ) · [ξ'' + U* · α' − a · α'']** (apparent mass group)
>       **− (2/μ) · [U*² · α + U* · ξ' + U* · (½ − a) · α']** (circulatory group)

> **Q_α = + (2/μ) · (½ + a) · [U*² · α + U* · ξ' + U* · (½ − a) · α']** (circulatory group)
>       **+ (1/μ) · [a · (ξ'' − a · α'') − U* · (½ − a) · α' − (1/8) · α'']** (apparent mass group)

Reading the structure:
- **U*² · α** terms → aerodynamic stiffness (∝ dynamic pressure)
- **U* · ξ'**, **U* · α'** terms → aerodynamic damping (linear in U*)
- **ξ''**, **α''** terms (no U*) → apparent-mass contributions

**note: incompressible derivation** when comp factor (PG) is added, its added into the 2π, so 2π becomes 2π/√(1−M∞²). Therefore, every circulatory term that born from 2π·α_eff **picks up an extra comp factor.** Hence, **circ term (quasi_steady.py) is multiplied by comp factor** Also dropped the apparent mass

Stage 1 is entirely circulatory terms.

### 3.4 Prandtl-Glauert compressibility correction

The lift-curve slope is scaled:

> 2π → 2π / √(1 − M_∞²)

This correction is applied *only to circulatory terms* (those originating from the 2π · α_eff steady contribution). The apparent-mass terms (those originating from π · (ξ'' + α' − a · α'') and the moment counterparts) are unaffected, as apparent mass is purely an inertial effect of the surrounding fluid.

In implementation, every coefficient 2 in front of a circulatory bracket in §3.3 becomes 2 / √(1 − M_∞²); coefficients in front of apparent-mass terms stay as written.

Validity: M_∞ < 0.7, strictly. Isogai Case A operates at M_∞ = 0.8, outside this range. PG is applied anyway as a baseline cheap correction, with the understanding that this is one of the inaccuracies the higher-fidelity stages will quantify.

### 3.5 Aero state contract

Quasi-steady has no memory: **n_aero_states = 0**. The `QuasiSteady` class implements:

- `n_aero_states = 0`
- `forces(tau, y_struct, y_aero, U_star)` returns [Q_ξ, Q_α] from §3.3 with PG correction from §3.4.
- `aero_rhs(...)` returns `np.zeros(0)`.

### 3.6 What this model captures and misses

**Captures correctly** (within thin-airfoil theory):
- Steady lift and moment.
- Static divergence (instability at U* = U*_D).
- Aerodynamic damping in plunge (positive damping below divergence).
- The qualitative phenomenology of sub-flutter decay → neutral stability → post-flutter growth → LCO (with β > 0).

**Does not capture**:
- Correct flutter speed. Quasi-steady lacks the wake-induced phase lag from C(k) that determines U*_F in classical flutter. Usually over-predicts U*_F for stiffness-dominated configurations.
- Correct frequencies of oscillation (errors growing with reduced frequency k). 
    - [reasonable at k <= 0.05, noticably off at k approx = 0.1, seriously wrong at k above 0.3]
- Transonic phenomena (shocks, the transonic dip in flutter speed). Prandtl-Glauert is a subsonic scaling, not a physics model.
- Stall and any viscous effects.

### 3.7 Steady-state sanity check

At ξ̇ = α̇ = ξ̈ = α̈ = 0 (frozen-in-time pose):

> C_L_QS → 2π · α / √(1 − M_∞²)
> C_M_QS → 2π · (½ + a) · α / √(1 − M_∞²)

This is the standard thin-airfoil, PG-corrected steady result. Any implementation must reduce to these at zero velocities. This is the first test the `QuasiSteady` class should pass.

### 3.8 Role in Multifidelity comparison

- cost accuracy floor: cheapest cost-accuracy model for aeroelastic studies
- interface validtion: first non-trivial implementation
- diagnostic baseline

---

## 4. Stage 2: Peters finite-state inflow

### 4.1 Assumption

Stage 1 set C(k) ≡ 1: the wake responds instantaneously, no memory, no phase lag. Stage 2 restores that memory. Theodorsen's C(k) is the *exact* frequency-domain wake lag, but it is transcendental (a ratio of Hankel functions) and cannot be written as a finite-state time-domain ODE. Peters, Karunamoorthy & Cao (1995) replace it with a **rational approximation**: the wake's induced velocity is expanded in N Glauert inflow states governed by a first-order ODE, whose transfer function (motion → induced velocity) approximates C(k). As N → ∞ the approximation converges to the exact C(k); N is an internal accuracy dial, with 3–8 states typically sufficient.

The induced-flow states are the *only* addition to the airfoil loading. Everything else — the circulatory lift slope, the apparent-mass terms — is unchanged in form from the underlying unsteady theory; Stage 2 adds the wake lag (and reinstates apparent mass, which Stage 1 dropped) on top of the same circulation.

### 4.2 Theory (Peters et al. 1995)

The N inflow states **Λ = [λ₁ … λ_N]ᵀ** evolve under a first-order ODE. In Peters' fully non-dimensional form (length on b, time on b/V), this is his Eq. (34):

> **[A]{λ̇} + {λ} = {c}(ẇ_a + ½ẇ_1)**

where the right-hand side is the airfoil's downwash forcing, and **[A]** is built from four N-only primitives (Eqs. 30, 31, 35, Appendix C):

> **D̄** (Eq. 30): bidiagonal, D̄_nm = +1/(2n) for n=m+1 (sub-diagonal), −1/(2n) for n=m−1 (super-diagonal), 0 otherwise.
> **c̄** (Eq. 31): c̄_n = 2/n.
> **d̄** (Eq. 31): d̄_n = ½ for n=1, else 0.
> **b̄** (Appendix C, augmented least squares): b̄_n = (−1)^(n−1)·(N+n−1)!/[(N−n−1)!·(n!)²] for n≠N; b̄_N = (−1)^(n−1).
> **[A] = D̄ + d̄b̄ᵀ + c̄d̄ᵀ + ½c̄b̄ᵀ** (Eq. 35).
    
D̄ is the bare state-to-state recursion; the three rank-one corrections are the algebraic price of the closure (below). The b̄ cap at n=N enforces Σb̄_n = 1, which cancels a 1/sinh η singularity in the closure (Eqs. C8–C9). Verification anchors: b̄ = [2,−1] (N=2), [6,−6,1] (N=3), [20,−90,140,−70,1] (N=5).

**Closure** (Eq. 28): the N states collapse to a single induced velocity at the airfoil,

> **λ₀ = ½ · (b̄ · Λ)**

### 4.3 Conversion to code convention

This codebase clocks time on ω_α (τ = ω_α·t), not on b/V. Under that re-scaling Peters' bare identity coefficient on {λ} becomes **U***, the reduced velocity (U* = U/(b·ω_α)) — Peters' hidden "1" is this code's U*, because the wake convection rate (V/b) reappears when time is normalised differently. The inflow ODE in code convention (`peters_finite.py: aero_rhs`) is:

> **A_bar · Λ' + U* · Λ = − c̄ · w'**,   with the 3/4-chord downwash rate  **w' = ξ'' + U* · α' + (½ − a) · α''**

The forcing w' is built entirely from structural motion (plunge acceleration, pitch rate, pitch acceleration) and is the abstract Glauert downwash of Peters' RHS specialised to the pitch–plunge DOFs.

**Effective angle of attack** gains exactly one term over Stage 1 — the induced velocity:

> **α_eff = α + ξ'/U* + (½ − a)·α'/U* + λ₀/U***

The circulatory bracket of §3.3 therefore gains a single term, `+ U*·λ₀`, with everything else (force prefactors, apparent-mass group) identical in form to Stage 1. Same forces, new angle of attack.

### 4.4 The EOM matrices

The aero generalized force is sorted by derivative order into the structural M/C/K matrices (`model/eom.py`, `model/analysis.py`). Two physical sources feed them — the circulatory bracket and the (reinstated) apparent-mass group:

- **K_aero** (∝ α), **C_aero** (∝ q̇): circulatory stiffness and damping — **identical to Stage 1 QS** (`peters_finite.py` reuses the QuasiSteady forms). First column of K_aero is zero (no aero force depends on plunge *position*).
- **M_a** (∝ q̈): apparent mass — `[[−1/μ, a/μ], [a/μ, −(a²+⅛)/μ]]`. Symmetric (genuine fluid inertia); carries no U* and no Prandtl–Glauert (inertia is speed/compressibility independent). Folded as `M_eff = M_s − M_a`. Reintroduced at Stage 2 (Stage 1 dropped it; a q̈-proportional aero force makes the EOM implicit, requiring it be folded into the mass matrix).
- **C_a** (∝ q̇): apparent damping — `[[0, U*/μ], [0, U*(½−a)/μ]]`. Only the pitch-rate (α') column is nonzero. Stored pre-flipped so it adds: `C_eff = C_s + C_aero + C_a`. The apparent-mass group contributes no stiffness (no displacement-proportional non-circulatory force), so there is no K_a.
- **K_aero_lambda** (∝ Λ): the inflow → structure coupling, `(2, N)`. It is the `U*·λ₀` slice of the circulatory bracket with the same row-prefactors, factored as `K_aero_lambda @ Λ`.
- **aero_forcing / aero_forcing_vel**: the structure → inflow coupling (RHS of the inflow ODE), split by derivative order — the q̈ slice `c̄ ⊗ [1, ½−a]` (E side) and the q̇ slice `c̄ ⊗ [0, −U*]` (A side).

These four blocks (structure↔structure, wake→structure, structure→wake, wake→wake) form the closed feedback loop that realises the wake lag: motion forces the wake (aero_forcing) → the wake evolves with a first-order lag (A_bar, U*) → its induced velocity λ₀ feeds force back to the structure (K_aero_lambda). Because Λ obeys a differential (not algebraic) equation, it can only chase the forcing at finite rate, arriving late — that lateness is the C(k) phase lag in the time domain.

### 4.5 Aero state contract

The state vector grows from 4 to **4 + N**: `[ξ, α, ξ', α', λ₁ … λ_N]`. The linear flutter analysis assembles a descriptor generalized eigenproblem `E·ẋ = A·x` (`analysis.py: descriptor_matrices`, `linearized_eigenvalues`); E is non-singular (det E = det(M_s − M_a)·det(A_bar) ≠ 0). Because the finite-state inflow makes the lift-deficiency function rational, the flutter analysis is a **p-method** (eigenvalues solved directly at each speed), not a classical k-method. The nonlinear LCO path time-marches the same 4+N system via `solver.py` (`solve_ivp`).

### 4.6 What this model captures and misses

**Captures:** the wake-lag phase that QS omits, hence a physical flutter boundary; apparent-mass inertia; arbitrary-motion response (the inflow ODE is not restricted to simple-harmonic motion).

**Misses:** still potential-flow, attached, thin-airfoil — no separation, no dynamic stall. No transonic physics (no shocks, no transonic dip, no Isogai Case A mechanism). 2D only. Finite-N truncation error (controllable via N). This is the boundary the Stage 4 CFD anchor exists to cross.

### 4.7 Validation

Validated against the Michigan pitch–plunge rig: predicted flutter speed **13.15 m/s vs 13.19 m/s experimental (0.3%)**, closing the −27% Stage 1 QS gap. Because the aerodynamic model is the only thing that changed between Stage 1 and Stage 2, this isolates the missing wake lag (not the structure) as the cause of the QS error. Correctness is further guarded by agreement between the two independent solver paths (eigenvalue/VGBF and nonlinear time-march), which agree on the flutter boundary and near-onset growth rate. LCO mechanics are validated (bounded, onset at the flutter speed, growth at the eigenvalue rate); quantitative LCO amplitude versus the García Pérez bifurcation diagram is a pending validation item.

### 4.8 Role in multifidelity comparison

Stage 2 is the cheapest model that gets the attached-flow flutter boundary physically right — it is the low-cost end of the cost-accuracy frontier for non-transonic regimes. Its single error source is **modelling error** (the rational C(k) approximation plus potential-flow assumptions); it has no discretisation error (it is analytical, with no mesh), which is the key structural contrast with the Stage 4 CFD model. N provides a cost-accuracy lever internal to Stage 2, before the much larger cost of the UVLM and CFD stages.

---

## 8. References

- Theodorsen, T. (1935). *General Theory of Aerodynamic Instability and the Mechanism of Flutter*. NACA Report 496.
- Fung, Y. C. (1955). *An Introduction to the Theory of Aeroelasticity*. Wiley.
- Bisplinghoff, R. L., Ashley, H., & Halfman, R. L. (1955). *Aeroelasticity*. Addison-Wesley.
- Lee, B. H. K., Price, S. J., & Wong, Y. S. (1999). Nonlinear aeroelastic analysis of airfoils: bifurcation and chaos. *Progress in Aerospace Sciences* 35, 205–334.
- Isogai, K. (1979). On the transonic dip mechanism of flutter of a sweptback wing. *AIAA Journal* 17(7), 793–795.
- Jones, R. T. (1940). *The Unsteady Lift of a Wing of Finite Aspect Ratio*. NACA Report 681.