## Eigenvalue analysis: from 2nd-order EOM to 1st-order state-space

### Why this section exists

The structural model's natural frequencies and damping ratios are computed via eigenvalue analysis. The naive approach — eigenvalues of M⁻¹K (a 2×2 problem) — works only in the special case of no damping and no aerodynamics. The general method requires converting the 2nd-order EOM into a 1st-order state-space system, then taking eigenvalues of the resulting 4×4 system matrix. This section documents why, how, and what the eigenvalues mean physically.

### The 2nd-order EOM and its limitation

The structural equation of motion is:

> M · q̈ + C · q̇ + K · q = 0

with q = [ξ, α]ᵀ, M, C, K all 2×2. This is a **2nd-order ODE in q**: the highest time derivative is q̈.

The naive eigenvalue analysis solves K·v = ω²·M·v (the generalized eigenvalue problem from the harmonic ansatz q ∝ e^(iωt) with C = 0). This returns 2 real eigenvalues ω² and 2 real eigenvectors. It works when:
- C = 0 (no damping)
- No aerodynamics (Q_aero = 0)
- Linearised around the trivial equilibrium (α_eq = 0)

For Stage 0 wind-off, all three conditions hold and the method gives ω₁ = 0.713, ω₂ = 5.34 (matching Isogai). But the method **cannot generalize** to damped or aerodynamic systems. The reasons are structural:

1. **C is not an input to `eigh(K, M)`** — there's nowhere to put it. Damping is invisible to the method.
2. **The 2×2 problem with real symmetric M and K is guaranteed real eigenvalues** (theorem of linear algebra). Real eigenvalues can only encode frequencies (oscillation rates), not damping (decay rates).
3. **Damping doesn't add new modes** — it shifts existing modes' eigenvalues from pure-imaginary to complex. The 2×2 method can't represent this shift.

### Why a 2-DOF system has 4 eigenvalues, not 2

For an N-DOF system, the state of the system at any instant is described by 2N variables (N positions + N velocities). The dynamics map this 2N-element state to its time derivative. Eigenvalue analysis on this map gives 2N eigenvalues.

For oscillatory systems, the 2N eigenvalues organize into N **complex-conjugate pairs**, one pair per physical mode. Each pair encodes one mode's frequency (Im λ) and damping rate (Re λ).

For N = 2 (your typical section): 4 eigenvalues = 2 pairs = 2 modes. **Two DOFs require four eigenvalues**, not two.

The 2×2 formula `K·v = ω²·M·v` returns only the *N* frequencies, throwing away the other *N* eigenvalue components (real parts) that encode damping. It works in the undamped case because the real parts are exactly zero — there's nothing to throw away. But this is a special case, not the general result.

### Converting to 1st-order state-space form

The conversion follows Wright & Cooper (2014), Eq. (10.16)–(10.20). The trick is to **augment the state vector** to include velocities as explicit state variables, and add the trivial kinematic identity that locks them to be the derivatives of position.

Define the augmented state:

> y = [ξ, α, ξ̇, α̇]ᵀ ∈ ℝ⁴

The trivial identity (Wright & Cooper Eq. 10.18) is:

> I · q̇ − I · q̇ = 0

This is vacuous as a statement of fact, but as a *matrix row* it serves a structural purpose: it gives the augmented system a row that tracks q̇ as a first-class state variable. Combined with the EOM, the 2nd-order system becomes 1st-order (Wright & Cooper Eq. 10.19–10.20):

> dy/dτ = Q · y

where Q is the **4×4 system matrix**:

> Q = [   0           I      ]
>     [ −M⁻¹·K     −M⁻¹·C    ]

Block structure: top-left and top-right are 2×2; bottom-left and bottom-right are 2×2. Total 4×4.

**What each block does:**

- **Top-right (I)**: encodes the kinematic identity "d/dτ of position = velocity." This is where q̇ enters the state explicitly.
- **Bottom-left (−M⁻¹·K)**: stiffness contribution to acceleration. From M·q̈ = −K·q solved for q̈.
- **Bottom-right (−M⁻¹·C)**: damping contribution to acceleration. From M·q̈ = −C·q̇ solved for q̈.
- **Top-left (0)**: position has no direct contribution to its own derivative — that's what velocity is for.

When aero is added (Stage 1+), M, C, K are augmented with aero contributions (M_aero, C_aero(U*), K_aero(U*)) before being assembled into Q. The block structure is unchanged.

### Why Q is non-symmetric (even when M, C, K are symmetric)

Each of M, C, K is real symmetric. But the block structure of Q is asymmetric: the top-right is I, while the bottom-left is −M⁻¹·K. These are different matrices, so Qᵀ ≠ Q.

Non-symmetric real matrices can have **complex eigenvalues** (no theorem forces them to be real, unlike the symmetric case). For oscillatory systems with damping, this is exactly what we want: complex eigenvalues encoding both frequency and damping rate.

The symmetry of M, C, K individually is preserved as a *physical* property (energy conservation in the undamped case, etc.), but the *eigenvalue problem* on Q is non-symmetric because of how the blocks are stacked. This non-symmetry is what unlocks the complex spectrum.

### Eigenvalues of Q: physical interpretation

For Q ∈ ℝ⁴ˣ⁴, the eigenvalue problem is:

> Q · v = λ · v

returning 4 complex eigenvalues λ₁, λ₂, λ₃, λ₄. For an oscillatory physical system, they organize as 2 complex-conjugate pairs:

> Mode 1: λ₁, λ₁* = σ₁ ± i·ω₁
> Mode 2: λ₂, λ₂* = σ₂ ± i·ω₂

(per Wright & Cooper Eq. 10.22). For each mode:

- **Frequency** = |Im λ| = ω
- **Damping rate** = −Re λ = −σ
- **Damping ratio** ζ = −Re λ / |λ|

Sign conventions:
- σ < 0 (Re λ < 0) → mode decays (stable)
- σ = 0 (Re λ = 0) → mode oscillates without decay (neutrally stable / flutter boundary)
- σ > 0 (Re λ > 0) → mode grows (unstable / post-flutter)

### Eigenvalue migration as a function of U*

For the wind-off undamped system (U* = 0, C = 0), all eigenvalues sit on the imaginary axis:

> λ = ±i·ω₁, ±i·ω₂

with real parts exactly zero. As U* increases (aerodynamics enters), the aero terms add damping and stiffness contributions to Q, and the eigenvalues **migrate** in the complex plane:

- Sub-flutter (U* < U*_F): aerodynamic damping pushes Re λ into the negative half-plane. Modes decay.
- At flutter (U* = U*_F): one mode's Re λ crosses zero. Sustained oscillation.
- Post-flutter (U* > U*_F): Re λ > 0 for that mode. Oscillation grows exponentially (linear theory); cubic stiffness β saturates it into LCO.

**Flutter prediction is then a single eigenvalue computation per U***. The 2×2 method cannot do this because it has no real-part axis to migrate along.

### Linearizing around α_eq ≠ 0

For the nonlinear stiffness model K(α), the eigenvalue analysis is done by *linearizing* K around an equilibrium pitch angle α_eq:

> K_linearized = K(α_eq) = K_α · (1 + β·α_eq²)

For α_eq = 0: standard linear analysis, gives Stage 0's ω₁ = 0.713, ω₂ = 5.34.
For α_eq ≠ 0: K_linearized has a different stiffness because of the cubic spring contribution. The eigenvalues shift accordingly.

This is one validation Steve requested: show that the frequencies vary with α_eq when β > 0. It validates both the nonlinear stiffness implementation *and* the eigenvalue machinery in one move.

### Relation to time-history integration

The state-space form `dy/dτ = Q·y` is the **linearized** form of the dynamics. The actual `structural_rhs` in `eom.py` computes the *nonlinear* dynamics (with state-dependent K via β·α²). Time-history integration via `solve_ivp` therefore captures the full nonlinear behavior, including amplitude-dependent frequency shifts.

The eigenvalue analysis trades nonlinearity for tractability: it answers "what does the dynamics look like in a neighborhood of α_eq?" In that neighborhood, the answer is exact (matrix eigenvalues are exact). Away from that neighborhood, the answer is approximate (the actual nonlinear system departs from the linearization).

Both views are useful: the linearized eigenvalues give you stability predictions at specific operating points (cheap, instant); the time-history gives you the full nonlinear behavior including post-flutter LCO (more expensive, but no linearization error).