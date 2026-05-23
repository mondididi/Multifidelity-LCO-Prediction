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

```
det(K - ω² M) = 0
det([1 - ω²     -1.8 ω²    ]) = 0
   ([-1.8 ω²    3.48 - 3.48 ω²])

(1 - ω²)(3.48 - 3.48 ω²) - (1.8 ω²)² = 0
3.48 (1 - ω²)² - 3.24 ω⁴ = 0
0.24 ω⁴ - 6.96 ω² + 3.48 = 0
```

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