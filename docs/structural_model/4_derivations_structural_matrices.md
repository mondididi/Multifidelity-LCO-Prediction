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