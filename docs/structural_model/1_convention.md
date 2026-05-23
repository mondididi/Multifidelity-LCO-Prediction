# 1. Convention


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

| Symbol  | Stage 0          | Stages 1–4                                                                 |
| ---     | ---              | ---                                                                        |
| Q_aero  | zero (no aero)   | computed by the aero model — possibly with extra internal aero states     |


## Sign conventions and gotchas

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
