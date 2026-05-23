## 2. Parameters

Eight nondimensional parameters, all stored on `TypicalSectionParams`:

| # | Symbol    | Code name     |IsogaiA| Meaning                                                           |
|---| ---       | ---           |---    | ---                                                               |
| 1 | a         | `a`           | −2.0  | Elastic-axis location, semichords from midchord, +aft             |
| 2 | x_α       | `x_alpha`     | 1.8   | Static imbalance — EA-to-CG distance, semichords, +aft            |
| 3 | r_α²      | `r_alpha_sq`  | 3.48  | Squared radius of gyration about the EA, normalised by semichord² |
| 4 | ω_h/ω_α   | `omega_ratio` | 1.0   | Uncoupled frequency ratio                                         |
| 5 | μ         | `mu`          | 60.0  | Mass ratio m / (π ρ b²)                                           |
| 6 | β         | `beta`        | 0.0   | Cubic stiffness coefficient in k_α(α) = k_α,0(1 + β α²)           |
| 7 | ζ_h       | `zeta_h`      | 0.0   | Plunge structural damping ratio                                   |
| 8 | ζ_α       | `zeta_alpha`  | 0.0   | Pitch structural damping ratio                                    |

Note: `a` is stored for completeness but doesn't appear in the structural
matrices — it's an aerodynamic quantity that enters the pitching-moment
calculation from Stage 1 onward.