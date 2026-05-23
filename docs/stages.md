# Stage Roadmap

This project compares aerodynamic models of increasing fidelity for predicting
nonlinear limit cycle oscillations (LCO) of a 2-DOF pitch–plunge typical
section. The structural model and validation benchmarks are held constant
across all stages; aerodynamic fidelity is the sole independent variable.


---


## Structural model (constant across all stages)

2-DOF pitch–plunge typical section with cubic pitch stiffness:

    k_α(α) = k_α,0 · (1 + β α²)

When β = 0 the system is linear; positive β bounds post-flutter growth into
LCO (Lee et al., *Prog. Aero. Sci.* 1999).

Canonical structural parameter set: Isogai (1979) NACA 64A010 Case A. See
`docs/structural_model/` for the full derivation and sanity checks.


---


## Stage 0 — Structural-only verification

**Goal:** Verify the structural model and time integrator in isolation, before
adding any aerodynamics.

**Aero loading:** none (`aero/none.py` returns zero force).

**Acceptance criteria:**
- Wind-off coupled natural frequencies of Isogai Case A match published
  values (0.713·ω_α and 5.34·ω_α) within 1%.
- Undamped free vibration with β ≠ 0 conserves total mechanical energy to
  relative drift < 10⁻⁶ over 20 nondimensional time units.

**Files:**
- `model/params.py`
- `model/eom.py`
- `model/solver.py`
- `model/analysis.py`
- `aero/base.py`, `aero/none.py`
- `validation/isogai_params.py`

**Deliverable:** `examples/stage0_free_vibration.py` — free-vibration time
history showing the two coupled modes.


---


## Stage 1 — Quasi-steady strip theory + Prandtl–Glauert

**Goal:** First aeroelastic predictions with the simplest possible aero model.

**Aero loading:** algebraic — instantaneous lift and pitching moment computed
from effective angle of attack. No memory, no aero states. Prandtl–Glauert
compressibility correction `1/√(1 − M∞²)` applied to C_Lα.

**Expected behaviour:** Will over-predict LCO amplitude and misplace the
bifurcation point relative to experiment. This is the *known weakness* the
higher-fidelity stages will improve on.

**File added:** `aero/quasi_steady.py`.

**Deliverable:** Flutter speed prediction vs Isogai's published linear flutter
boundary; LCO bifurcation diagram (amplitude vs U*) for β = 3.


---


## Stage 2 — Wagner indicial + Jones two-pole approximation

**Goal:** Capture the unsteady delay in circulatory lift build-up.

**Aero loading:** Wagner's indicial response approximated by Jones' two-pole
rational form. Introduces two aerodynamic state variables; the integrator
state grows from 4 to 6 components.

**Expected behaviour:** Significant improvement in LCO amplitude prediction
over quasi-steady, at modest computational cost.

**File added:** `aero/wagner.py`.

**Deliverable:** Comparison of LCO bifurcation diagrams Stage 1 vs Stage 2 vs
experiment.


---


## Stage 3 — Unsteady Vortex Lattice Method (UVLM)

**Goal:** Capture wake roll-up and near-field vortex interactions without the
expense of solving the Navier–Stokes equations.

**Aero loading:** discrete vortices shed into a wake mesh at each timestep.
Aerodynamic state is the vector of wake-vortex circulations (potentially
hundreds of states).

**Expected behaviour:** Physically richer than Wagner; required to capture
wake-airfoil interactions. Captures attached-flow LCO well; cannot resolve
shocks or boundary-layer separation.

**File added:** `aero/uvlm.py`.

**Deliverable:** UVLM bifurcation diagram; comparison of all three low-order
stages against experiment.


---


## Stage 4 — CFD coupling

**Goal:** Reach the regime where CFD is genuinely required (transonic dip,
nonlinear shock motion).

Two implementations to enable an intermediate-vs-high CFD cost comparison:

1. **Intermediate-cost method** (choice pending — discuss with Abhijith):
   harmonic balance, linearised frequency-domain Euler, or POD/HDHB
   reduced-order model.
2. **High-cost method:** unsteady Euler with RBF mesh morphing (Yuan, Sandhu
   & Poirel, *J. Aerospace Eng.* 2021).

**Expected behaviour:** Both methods resolve the transonic dip near Isogai's
design point; the intermediate-cost method should approach the unsteady-Euler
result at a fraction of the cost.

**Files added:** `aero/cfd_intermediate.py`, `aero/cfd_euler.py`.

**Deliverable:** Cost-accuracy frontier map across all five stages, organised
by flow regime (attached / transonic / stall-bounded).


---


## Validation targets

| Target | Source | Used for |
|---|---|---|
| Isogai (1979) linear flutter boundary             | Theoretical   | Stage 0 eigen-check; Stage 1+ linear-flutter cross-check                                                                  |
| Lee et al. (1999) post-Hopf LCO data              | Numerical     | Cubic-spring LCO amplitude/frequency for low-order stages                                                                 |
| AGARD CT6 / Davis (1982)                          | Experimental  | CFD prescribed-motion validation before aeroelastic coupling (Stage 4)                                                    |
| Bristol LCO rig — Tartaruga et al. (IFASD 2019)   | Experimental  | Closed-loop LCO validation; *note: freeplay nonlinearity, not cubic*                                                      |
| García Pérez et al. (2024), Michigan/Toulouse     | Experimental  | Closed-loop LCO validation with cubic-spring nonlinearity — exact structural match. To confirm with Fintan whether to add |


---


## Constants held fixed across all stages

To isolate aerodynamic fidelity as the sole independent variable:

- Airfoil geometry: NACA 64A010
- Structural parameters: Isogai Case A
- Cubic stiffness coefficient β
- Time-step size and integrator tolerances
- Validation targets

Any change to the above mid-project invalidates earlier comparisons. Lock these
before Stage 1.