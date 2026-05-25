# Structural Model Reference

A standalone reference for the 2-DOF pitch–plunge structural model implemented
in `mflco.model`. Covers the equations of motion, the structural matrices, the
nondimensionalisation, and why each matrix has the form it does.

If anything in `params.py`, `eom.py`, or `analysis.py` ever looks wrong, start
here.


## 1. Convention

**State vector** (nondimensional time τ = ω_α · t, with ω_α the uncoupled pitch
natural frequency):

    y = [ξ, α, ξ', α']ᵀ

| Symbol     | Meaning                                                            |
| ---        | ---                                                                |
| ξ = h/b    | Plunge displacement, normalised by semichord b. Positive downward. |
| α          | Pitch angle, radians. Positive nose-up.                            |
| ξ' = dξ/dτ | Nondimensional plunge velocity.                                    |
| α' = dα/dτ | Nondimensional pitch rate.                                         |

**Governing equation** (matrix form):

    M q̈ + C q̇ + K(q) q = Q_aero(τ, y)

where q = [ξ, α]ᵀ and dots denote derivatives with respect to τ. Stage 0 sets
Q_aero = 0; later stages provide it from the aero model.

**Definitions of terms in the governing equation:**

| Symbol     | Type             | Meaning                                                                                                                                            |
| ---        | ---              | ---                                                                                                                                                |
| q          | 2-vector         | Displacement vector `[ξ, α]ᵀ`. The two structural DOFs.                                                                                            |
| q̇          | 2-vector         | Velocity vector `[ξ', α']ᵀ`. First time-derivative of q.                                                                                           |
| q̈          | 2-vector         | Acceleration vector `[ξ'', α'']ᵀ`. Second time-derivative of q.                                                                                    |
| M          | 2×2 matrix       | Structural mass matrix. Symmetric, constant. Encodes inertia and the EA–CG offset coupling. See §4.1.                                              |
| C          | 2×2 matrix       | Structural damping matrix. Diagonal, constant. Encodes plunge and pitch damping independently. See §4.2.                                           |
| K(q)       | 2×2 matrix       | Structural stiffness matrix. Diagonal, **α-dependent** through the cubic spring `k_α(1 + β α²)`. Reduces to constant when β = 0. See §4.3.         |
| Q_aero     | 2-vector         | Generalised aerodynamic force vector `[Q_ξ, Q_α]ᵀ`. Nondimensional lift and pitching-moment contributions to the plunge and pitch DOFs respectively. |
| τ          | scalar           | Nondimensional time, `τ = ω_α · t`. Time scaled by the uncoupled pitch natural frequency.                                                          |
| y          | 4-vector         | Full state vector `[ξ, α, ξ', α']ᵀ`. What the integrator advances in time.                                                                         |

**Stage-dependent values:**

| Symbol  | Stage 0          | Stages 1–4                                                             |
| ---     | ---              | ---                                                                    |
| Q_aero  | zero (no aero)   | computed by the aero model — possibly with extra internal aero states  |


### Sign conventions and gotchas

- **Pitch sign:** α > 0 means nose-up. The cubic spring restoring force is
  always opposite to α (the spring resists deflection in either direction).
  Make sure β never goes negative in code unless you're deliberately modelling
  a softening spring (rare in this project).
- **Plunge sign:** h > 0 (and ξ > 0) means downward by the convention used in
  Lee et al. and most LCO literature. If you compare against a paper using
  upward-positive plunge, the lift-force sign flips.
- **Elastic-axis location a:** a = −2.0 means the EA is 2 semichords *upstream*
  of midchord — i.e. forward of the leading edge. Sounds non-physical for a 2D
  section, but it's deliberate (mimics outboard section of a swept-back wing).
  Don't "fix" it.
- **Mass-matrix off-diagonals:** must be equal and equal to x_α. If you ever
  see `[[1, x_α/r_α²], [x_α, 1]]` somewhere, that's the asymmetric form leaking
  in. Convert it back.


## 2. Parameters

Eight nondimensional parameters stored on `TypicalSectionParameters`. Defaults
reproduce Isogai (1979) Case A — a swept-back wing outboard section at
transonic conditions.

| # | Symbol  | Code name     | Isogai A | Meaning                                                           |
|---| ---     | ---           | ---      | ---                                                               |
| 1 | a       | `a`           | −2.0     | Elastic-axis location, semichords from midchord, +aft             |
| 2 | x_α     | `x_alpha`     | +1.8     | Static imbalance — EA-to-CG distance, semichords, +aft            |
| 3 | r_α²    | `r_alpha_sq`  | 3.48     | Squared radius of gyration about the EA, normalised by b²         |
| 4 | ω_h/ω_α | `omega_ratio` | 1.0      | Uncoupled frequency ratio (plunge over pitch)                     |
| 5 | μ       | `mu`          | 60.0     | Mass ratio m / (π ρ b²)                                           |
| 6 | β       | `beta`        | 0.0      | Cubic stiffness coefficient in k_α(α) = k_α,0(1 + β α²)           |
| 7 | ζ_h     | `zeta_h`      | 0.0      | Plunge structural damping ratio                                   |
| 8 | ζ_α     | `zeta_alpha`  | 0.0      | Pitch structural damping ratio                                    |


### Per-parameter notes

- **a** is stored for completeness but doesn't appear in the structural
  matrices — it's an aerodynamic quantity that enters the pitching-moment
  calculation from Stage 1 onward.
- **x_α** drives the inertia coupling between plunge and pitch via the
  off-diagonal of M. System decouples completely when x_α = 0 (see §6.1).
- **ω_h/ω_α = 1.0** in Isogai Case A is deliberate worst-case — equal
  uncoupled frequencies make the two modes easiest to couple. After full
  coupling, wind-off frequencies split into 0.713·ω_α and 5.34·ω_α (§6.2).
- **μ = 60** is a typical aircraft-wing value. Higher μ → heavier structure
  relative to air → slower flutter onset. Lower μ → stronger aero coupling.
- **β** is zero by default, recovering Isogai's linear case (used for Stage 0
  verification). β > 0 produces a hardening pitch spring that bounds
  post-flutter growth into a limit cycle. Typical experimental values β ≈ 3
  to 100 (Lee et al. 1999, §5.2).
- **ζ_h, ζ_α** are zero by default following Isogai. Real wings have ζ ≈
  0.01–0.05; turn on for comparisons against experimental data.


## 3. Why M is symmetric and **not** diagonal

This is the most common point of confusion. The mass matrix has off-diagonal
terms; that is correct and necessary.


### 3.1 Physical origin

The elastic axis (EA, where the springs act) and the centre of mass (CG, where
inertia is located) are **not at the same point**. They are separated by a
distance `x_α · b`. For Isogai Case A, x_α = 1.8 — the CG sits 1.8 semichords
aft of the EA.

Because of this offset:
- A *pitch acceleration* α̈ produces a *linear* acceleration of the CG, which
  shows up as a force in the plunge equation.
- A *plunge acceleration* ḧ produces a *torque* about the EA via the offset CG,
  which shows up as a moment in the pitch equation.

The two DOFs are **mechanically coupled through inertia**. The off-diagonal of M
encodes that coupling.


### 3.2 Dimensional derivation (Lee et al. 1999, Eqs. 8–9)

Newton's second law about the EA, dimensional form:

    m ḧ + S_α α̈ + K_h h = −L
    S_α ḧ + I_α α̈ + K_α α = M_ea

where:
- m = section mass per unit span
- S_α = m · x_α · b = static mass moment about the EA (this is the off-diagonal coupling)
- I_α = moment of inertia about the EA = m · r_α² · b²
- K_h, K_α = linear plunge and pitch spring stiffnesses
- L, M_ea = aerodynamic lift and moment about the EA

In matrix form:

    [m    S_α] [ḧ]    [K_h  0 ] [h]   [-L  ]
    [S_α  I_α] [α̈]  + [0    K_α] [α] = [M_ea]

Mass matrix is symmetric (both off-diagonals are S_α) but **not** diagonal
unless S_α = 0, i.e. unless the CG coincides with the EA.


### 3.3 Nondimensionalising — the first attempt gives an asymmetric M

Scale time by ω_α (so τ = ω_α t), and the plunge by the semichord (ξ = h/b).
Divide the first equation by `m b ω_α²` and the second by `I_α ω_α²`:

    First eq:  ξ̈ + x_α α̈ + (ω_h/ω_α)² ξ = -L / (m b ω_α²)
    Second eq: (x_α / r_α²) ξ̈ + α̈ + α = M_ea / (I_α ω_α²)

Using the identities S_α / (m b) = x_α, K_h/(m ω_α²) = (ω_h/ω_α)², and
K_α / (I_α ω_α²) = 1. Notice the mass matrix now reads:

    M_asym = [1          x_α]
             [x_α/r_α²    1 ]

Symmetric? **No** — the (0,1) entry is `x_α` but the (1,0) entry is `x_α/r_α²`.
That's the price of dividing the two equations by different normalising factors.


### 3.4 The symmetric Lagrangian form

To recover symmetry, multiply the *second equation* through by `r_α²`. This
scales the entire pitch row uniformly — it doesn't change the dynamics, just
rescales how that equation is written. The result:

    M = [1     x_α  ]
        [x_α   r_α² ]

Symmetric. This is what `mass_matrix()` returns in `params.py`. Equivalent
multiplication is applied to the damping and stiffness matrices:

    C = [2 ζ_h (ω_h/ω_α)    0          ]
        [0                  2 ζ_α r_α² ]

    K(α) = [(ω_h/ω_α)²    0                  ]
           [0             r_α² (1 + β α²)    ]

Note that K[1,1] now carries an `r_α²` factor that it wouldn't have in the
asymmetric form — this is purely the consequence of the row-rescaling. **The
eigenfrequencies of (M, K) are identical to those of the asymmetric form**,
because multiplying one row by a constant doesn't change the generalised
eigenvalues. Same physics, cleaner form.


### 3.5 Why bother with symmetry?

Two reasons:

1. **Energy is well-defined.** For a Lagrangian system, kinetic energy is
   KE = ½ q̇ᵀ M q̇. This formula only gives a physically correct (positive
   definite, conserved) energy when M is symmetric. Without symmetry, the same
   formula gives a number that drifts during integration even with a perfect
   integrator — making energy conservation useless as a verification tool.

2. **Eigenproblems are cleaner.** Generalised eigenproblems `K x = ω² M x` with
   symmetric M and K have real eigenvalues and orthogonal eigenvectors. The
   asymmetric form has the same eigenvalues but the algebra is messier.

The energy-conservation test in `test_eom_nonlinear.py` requires the symmetric
form. If you ever see that test failing with large drift, it usually means
someone has accidentally reverted M to the asymmetric form.


## 4. The three structural matrices

This section gives the final matrix forms and explains where each entry comes
from. The full derivation is in §4.0; the three subsections that follow (§4.1,
§4.2, §4.3) summarise the resulting M, C, K with their per-entry notes.


### 4.0 Derivation: dimensional → asymmetric nondimensional → symmetric

All three matrices follow the same derivation. Damping and stiffness behave
analogously to mass; the only reason they look different is the row-rescaling
in step 3.

**Step 0 — Dimensional EOM** (Dowell §3.2.1, p. 63, Eq. 3.2.15, with damping
and cubic spring added):

    m ḧ   + S_α α̈ + C_h ḣ + K_h h             = -L
    S_α ḧ + I_α α̈ + C_α α̇ + K_α (α + β α³)   = M_ea

In matrix form:

    [m    S_α] [ḧ]    [C_h  0 ] [ḣ]    [K_h  0            ] [h]   [-L  ]
    [S_α  I_α] [α̈] +  [0    C_α] [α̇] + [0    K_α(1 + β α²)] [α] = [M_ea]

All three dimensional matrices are already symmetric. M has off-diagonal
coupling through S_α; C and K are diagonal.

**Step 1 — Introduce damping ratios** (Lee et al. 1999, p. 223, definitions
above Eq. 10; standard `2ζω` convention):

    C_h = 2 m ζ_h ω_h           ζ_h = C_h / (2 m ω_h)
    C_α = 2 I_α ζ_α ω_α         ζ_α = C_α / (2 I_α ω_α)

where ω_h = √(K_h/m) and ω_α = √(K_α/I_α) are the uncoupled natural frequencies.

**Step 2 — Nondimensionalise** with ξ = h/b and τ = ω_α t. Substitute time
derivatives (ḣ = b ω_α ξ', ḧ = b ω_α² ξ'', etc.), divide the first equation
by `m b ω_α²` and the second by `I_α ω_α²`. Using the identities

    S_α / (m b) = x_α           K_h / (m ω_α²) = (ω_h/ω_α)²
    S_α b / I_α = x_α / r_α²    K_α / (I_α ω_α²) = 1

gives the asymmetric nondimensional form:

    Plunge:  ξ''            + x_α α''  + 2 ζ_h (ω_h/ω_α) ξ' + (ω_h/ω_α)² ξ   = Q_ξ
    Pitch:   (x_α/r_α²) ξ'' + α''      + 2 ζ_α α'           + (1 + β α²) α   = Q_α

Reading off the matrices:

    M_asym = [1         x_α]   C_asym = [2 ζ_h (ω_h/ω_α)  0     ]
             [x_α/r_α²   1 ]            [0                 2 ζ_α ]

    K_asym(α) = [(ω_h/ω_α)²   0          ]
                [0             1 + β α²  ]

Mass matrix is asymmetric: M[0,1] = x_α but M[1,0] = x_α/r_α². Damping and
stiffness are diagonal — the asymmetry only affects M at this stage.

**Step 3 — Symmetrise by multiplying the pitch equation by r_α²**.
This is a single operation applied to *every* term in the second equation:

    (x_α/r_α²) ξ''   · r_α²   →   x_α ξ''
    α''              · r_α²   →   r_α² α''
    2 ζ_α α'         · r_α²   →   2 ζ_α r_α² α'
    (1 + β α²) α     · r_α²   →   r_α² (1 + β α²) α
    Q_α              · r_α²   →   r_α² Q_α

Mathematically equivalent (multiplying both sides by a nonzero constant
preserves solutions), but now M is symmetric. Every entry in row 2 of all
three matrices picks up a factor of r_α²:

| Position | Before (asym) | After symmetrise   |
| ---      | ---           | ---                |
| M[1,0]   | x_α / r_α²    | **x_α**            |
| M[1,1]   | 1             | **r_α²**           |
| C[1,1]   | 2 ζ_α         | **2 ζ_α r_α²**     |
| K[1,1]   | 1 + β α²      | **r_α² (1 + β α²)**|

The first equation (plunge) is unchanged.

**Three patterns worth remembering:**

1. `ω_h/ω_α` appears in plunge entries because plunge oscillates at ω_h but
   we scale time by ω_α. The ratio is the nondimensional plunge frequency.
2. `r_α²` appears in pitch entries because we multiplied the pitch equation
   by r_α² to symmetrise M. Anything in row 2 inherits this factor.
3. `(1 + β α²)` is the spring law itself, separate from any scaling. It would
   be there in dimensional form too — it's the physics of the cubic spring.


### 4.1 Mass matrix

    M = [1     x_α  ]
        [x_α   r_α² ]

- **Symmetric**, not diagonal. Off-diagonal = x_α (the inertial coupling).
  Becomes diagonal only in the trivial limit x_α = 0, which decouples plunge
  and pitch entirely (see §6.1).
- Independent of α — mass doesn't change with motion.
- Returned by `params.mass_matrix()`.
- Constant throughout a simulation; could be cached if performance ever matters.


### 4.2 Damping matrix
(Standard 2ζω convention; Lee et al. 1999 p. 223 damping-ratio definitions
above Eq. (10); row-2 r_α² factor from the row-rescaling in §4.0 step 3.)

    C = [2 ζ_h (ω_h/ω_α)   0          ]   ← plunge equation (cf. Lee Eq. 10)
        [0                 2 ζ_α r_α² ]   ← pitch equation  (cf. Lee Eq. 11, ×r_α²)

- **Diagonal.** Plunge and pitch dampers are independent in this model.
- Default values ζ_h = ζ_α = 0 → C is the zero matrix.
- Each entry follows the standard `2ζω` form from textbook structural dynamics.
  For plunge, ω is the nondimensional plunge frequency ω_h/ω_α. For pitch, ω
  is ω_α/ω_α = 1, so only `2 ζ_α` would appear if not for the r_α² factor
  from the row rescaling.
- Returned by `params.damping_matrix()`.


### 4.3 Stiffness matrix
(Cubic spring form: proposal Eq. (1); Lee et al. 1999 §5.2, p. 256; row-2
r_α² factor from the row-rescaling in §4.0 step 3.)

    K(α) = [(ω_h/ω_α)²   0                ]
           [0            r_α² (1 + β α²)  ]

- **Diagonal**, but **α-dependent** through the cubic spring `k_α(1 + β α²)`.
- K[0,0] = (ω_h/ω_α)² — linear plunge spring, from K_h / (m ω_α²).
- K[1,1] = r_α² (1 + β α²) — pitch spring including hardening. The leading
  r_α² comes from the row rescaling that gives M its symmetric form.
- When β = 0, K[1,1] = r_α² and the system is linear; this is Isogai's
  original case.
- When β > 0, the spring hardens with amplitude — this is the mechanism that
  bounds post-flutter growth into a limit cycle instead of unbounded divergence
  (see Lee et al., Prog. Aero. Sci. 1999, §5.2 p. 256).
- Returned by `params.stiffness_matrix(alpha)`. Note the argument — must pass
  the current α because the matrix depends on it.


## 5. Energy of the system

For the symmetric Lagrangian form:

**Kinetic energy:**

    KE = ½ q̇ᵀ M q̇
       = ½ (ξ'² + 2 x_α ξ' α' + r_α² α'²)

**Potential energy** (linear plunge spring + nonlinear pitch spring):

    PE = ½ (ω_h/ω_α)² ξ²              (plunge)
       + ½ r_α² α² + ¼ r_α² β α⁴      (pitch, including cubic term)

The cubic-spring potential comes from integrating the pitch restoring force:

    ∫ k_α,0 (1 + β α²) α dα = k_α,0 (½ α² + ¼ β α⁴)

and the leading `r_α²` factor again comes from the row rescaling. Note: the
nondimensionalisation absorbs `k_α,0` into the unit choice (K_α/(I_α ω_α²) = 1).

**Total energy** E = KE + PE is conserved for undamped (ζ = 0) free vibration
(Q_aero = 0). The test in `test_eom_nonlinear.py` verifies this conservation to
relative drift < 10⁻⁶ over 20 nondimensional time units.


## 6. Sanity checks

Quick mental checks against which to validate any reimplementation.


### 6.1 Uncoupled limit: x_α = 0

Set x_α = 0 (CG coincident with EA). Then:

- M = diag(1, r_α²) — diagonal.
- K = diag((ω_h/ω_α)², r_α²) — already diagonal.
- Eigenproblem: ω₁² = (ω_h/ω_α)², ω₂² = 1.
- So the system reduces to two independent oscillators at the uncoupled
  natural frequencies. No flutter possible, no LCO possible.

This is tested in `test_uncoupled_limit_recovers_input_frequencies`.


### 6.2 Isogai Case A: published wind-off frequencies

For the Isogai parameters (x_α = 1.8, r_α² = 3.48, ω_h/ω_α = 1, β = 0):

    det(K - ω² M) = 0
    det([1 - ω²     -1.8 ω²        ]) = 0
       ([-1.8 ω²    3.48 - 3.48 ω²])

    (1 - ω²)(3.48 - 3.48 ω²) - (1.8 ω²)² = 0
    3.48 (1 - ω²)² - 3.24 ω⁴ = 0
    0.24 ω⁴ - 6.96 ω² + 3.48 = 0

Solving:
- ω₁² ≈ 0.5089 → ω₁ ≈ 0.7134
- ω₂² ≈ 28.49 → ω₂ ≈ 5.338

These match Isogai's published values 0.713 and 5.34 to three significant
figures. Tested in `test_isogai_case_a_eigenfrequencies`.


### 6.3 Linear limit: β = 0

With β = 0, K[1,1] reduces to r_α² regardless of α. The stiffness matrix
becomes independent of α, and the system is linear. Tested in
`test_stiffness_matrix_is_linear_when_beta_zero`.


### 6.4 Hardening: β > 0

For β > 0, K[1,1] grows with α². At α = 0.3 rad with β = 3, the multiplier
becomes 1 + 3 × 0.09 = 1.27. Tested in
`test_stiffness_matrix_hardens_when_beta_positive`.


### 6.5 Energy conservation

Undamped (ζ = 0), no aero, β ≠ 0, small pitch perturbation. Run the integrator
for 20 nondimensional time units with tight tolerances. Total energy E should
not drift by more than 1 part in 10⁶. Tested in
`test_energy_conservation_undamped_nonlinear`.


## 7. References

- Dowell, E. H. (ed.) (2015). *A Modern Course in Aeroelasticity*, 5th ed.,
  Springer. §3.2.1 pp. 61–64 for the Lagrangian derivation of the typical
  section EOM.
- Isogai, K. (1979). "On the transonic-dip mechanism of flutter of a swept-back
  wing." *AIAA Journal* 17(7), 793–795. Original Case A parameter set.
- Lee, B.H.K., Price, S.J., Wong, Y.S. (1999). "Nonlinear aeroelastic analysis
  of airfoils: bifurcation and chaos." *Progress in Aerospace Sciences* 35,
  205–334. Authoritative treatment of the cubic-spring LCO formulation;
  post-Hopf bifurcation data for various μ, x_α, r_α², β.
- Yuan, Y., Sandhu, R., Poirel, D. (2021). "Fully Coupled CFD/CSD Aeroelastic
  Simulations of a Flexible Wing with Underwing Stores." *J. Aerospace Eng.*
  34(2): 04020117. Cleanly reproduces Isogai's Case A parameters and the
  wind-off coupled frequencies (0.713 ω_α, 5.34 ω_α).
- García Pérez, J., Ghadami, A., Sanches, L., Epureanu, B. I., Michon, G.
  (2024). "Data-Driven Bifurcation Analysis of Experimental Aeroelastic
  Systems Using Preflutter Measurements." *AIAA Journal* 62(5), 1906–1915.
  Experimental 2-DOF airfoil with cubic spring; potential validation target
  alongside or instead of Bristol/Tartaruga (2019).