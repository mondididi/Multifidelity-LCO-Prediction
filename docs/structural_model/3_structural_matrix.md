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


### 3.2 Dimensional derivation: **Lee(1999) (8,9) (ref)

Newton's second law about the EA, dimensional form:

```
m ḧ + S_α α̈ + K_h h = −L
S_α ḧ + I_α α̈ + K_α α = M_ea
```

where:
- m = section mass per unit span
- S_α = m · x_α · b = static mass moment about the EA (this is the off-diagonal coupling)
- I_α = moment of inertia about the EA = m · r_α² · b²
- K_h, K_α = linear plunge and pitch spring stiffnesses
- L, M_ea = aerodynamic lift and moment about the EA

In matrix form:

```
[m    S_α] [ḧ]  +  [K_h  0  ] [h]   [-L  ]
[S_α  I_α] [α̈]  +  [0    K_α] [α] = [M_ea]
```

Mass matrix is symmetric (both off-diagonals are S_α) but **not** diagonal
unless S_α = 0, i.e. unless the CG coincides with the EA.


### 3.3 Nondimensionalising — the first attempt gives an asymmetric M

Scale time by ω_α (so τ = ω_α x t), and the plunge by the semichord (ξ = h/b).
Divide the first equation by `m b ω_α²` and the second by `I_α ω_α²`:

```
First eq:  ξ̈ + x_α α̈ + (ω_h/ω_α)² ξ = -L / (m b ω_α²)
Second eq: (x_α / r_α²) ξ̈ + α̈ + α = M_ea / (I_α ω_α²)
```

Using the identities S_α / (m b) = x_α, K_h/(m ω_α²) = (ω_h/ω_α)², and
K_α / (I_α ω_α²) = 1. Notice the mass matrix now reads:

```
M_asym = [1          x_α]
         [x_α/r_α²    1 ]
```

Symmetric? **No** — the (0,1) entry is `x_α` but the (1,0) entry is `x_α/r_α²`.
That's the price of dividing the two equations by different normalising factors.


### 3.4 The symmetric Lagrangian form

To recover symmetry, multiply the *second equation* through by `r_α²`. This
scales the entire pitch row uniformly — it doesn't change the dynamics, just
rescales how that equation is written. The result:

```
M = [1     x_α  ]
    [x_α   r_α² ]
```

Symmetric. This is what `mass_matrix()` returns in `params.py`. Equivalent
multiplication is applied to the damping and stiffness matrices:

```
C = [2 ζ_h (ω_h/ω_α)        0          ]
    [0                       2 ζ_α r_α² ]

K(α) = [(ω_h/ω_α)²    0                    ]
       [0             r_α² (1 + β α²)     ]
```

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