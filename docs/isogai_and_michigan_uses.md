# Validation strategy: why both Isogai and Michigan

## The distinction in one line

- **Isogai Case A — *verification*.** Answers "**is my code correct?**" against a known, published, exactly-defined benchmark.
- **García Pérez et al. (Michigan) — *validation*.** Answers "**does my model match reality?**" against a physical wind-tunnel experiment.

These are two different questions. The project needs both, and they cannot be collapsed into one.

## Quick reference: which reference, and when

| Reference                 | Type                              | Role                                              | When it is used                                                               |
|---                        |---                                |---                                                |---                                                                            |
| Isogai Case A             | Numerical benchmark               | **Verification** — "is the code correct?"         | At *every* fidelity rung, before that rung is trusted on real data            |
| García Pérez (Michigan)   | Physical experiment               | **Validation** — "does the model match reality?"  | After a rung is verified; the headline quantitative experimental comparison   |
| Bristol Tartaruga         | Physical experiment (secondary)   | **Qualitative cross-check only**                  | Optional sanity that behaviour resembles another real pitch-plunge rig        |

Three rules that go with the table:

- **Verification repeats per rung.** Isogai does not verify the skeleton once and for all. Each aero model (quasi-steady, Peter's, UVLM, CFD) gets its own Isogai pass before it is compared to Michigan. "Skeleton verified" is not "QS aero verified" — they are separate confirmations.
- **Calibrate the structure, never the aero.** When building the Michigan section, calibrate only the *structural / sectional* properties to the rig's *measured modal data* (5.3 / 6.2 Hz, structural damping). The published components do not self-consistently reproduce the measured modes, so tuning the effective inertia and stiffness ratios is legitimate model setup. Never tune the *aero* model to make flutter speed or LCO amplitude match the experiment — the aero stays fixed, and the gap between its prediction and reality *is the result*. Fitting the aero to the answer destroys the cost-accuracy comparison.
- **Tartaruga is qualitative only.** The Bristol rig does not report frequency, damping, and LCO-amplitude data in a form that can be overlaid quantitatively, so it serves at most as a "does our qualitative behaviour look like a real pitch-plunge LCO" cross-check. Michigan is the quantitative experimental anchor.

## Why Isogai earns its place

Isogai Case A is not a physical rig — it is a *defined* nondimensional benchmark (`a = −2.0`, `x_alpha = 1.8`, `r_alpha_sq = 3.48`, `mu = 60`, `omega_ratio = 1.0`). Three consequences follow:

1. **No parameter extraction.** The numbers *are* the section. There is no digitising, no unit conversion, no unmodelled mass — so nothing upstream of the solver can be wrong.
2. **A known answer.** Its flutter behaviour is published and replicated across dozens of papers. If the solver disagrees with Isogai, the cause is *unambiguously the solver*. It is a unit test with the answer printed in the back of the book.
3. **Full eigenvalue data and zero structural damping.** Frequencies *and* dampings are available, and the structure is undamped by design, so the flutter boundary is a pure property of the aerodynamics. A discrepancy points straight at the aero implementation, not at a damping assumption.

**Heads-up 2** — what "verify on Isogai" means at the QS rung. This is a CONSISTENCY
check (U*=0 oracle passes, modes behave sensibly, the expected QS over-
destabilisation appears). It is NOT "match a published Isogai flutter speed":
Isogai's famous result is the transonic-dip boundary computed at high fidelity
(Euler/CFD), which is the target for the CFD rung -- not for QS. QS-on-Isogai
verifies the plumbing and demonstrates the known QS limitation; it is not
expected to hit a specific Isogai flutter number.

## Why Michigan cannot do the verification job

García Pérez is a real rig, which is exactly what makes it the right *validation* target — and the wrong *verification* tool:

- **Confounded failure modes.** If a Michigan prediction is wrong, the cause could be (a) a code bug, (b) a parameter-extraction error, or (c) unmodelled physics in the rig. All three are entangled, so you can never cleanly attribute a disagreement.
- **No damping data.** Figure 6 reports natural frequencies only — there is no experimental damping curve, so the damping branch of a VGBF cannot be validated against it.
- **The rig is not even self-consistent on paper.** Assembling the published component values (m = 2.4 kg, I_α = 3.401×10⁻³ kg·m², K_h = 1580 N/m, K_α = 11.53 N·m/rad) gives coupled wind-off frequencies of **~4.0 and 9.9 Hz**, not the **measured 5.3 and 6.2 Hz**. The effective pitch inertia is roughly twice the bare-wing value once the spar and spring assembly are included, so the nondimensional parameters must be *calibrated* to the measured modes rather than computed bottom-up. If Michigan were the only reference, a wrong U\*=0 frequency could not be distinguished from a code bug.

## Division of labour

## Division of labour

|                       | Isogai Case A                         | García Pérez (Michigan)               |
|---                    |---                                    |---                                    |
| Role                  | Verification                          | Validation                            |
| Question answered     | Is the code correct?                  | Does the model match reality?         |
| Source of truth       | Published nondimensional benchmark    | Physical wind-tunnel rig              |
| Parameters            | Defined exactly                       | Extracted / calibrated from the rig   |
| Structural damping    | Zero (aero isolated)                  | Real, must be estimated               |
| Frequencies           | Known                                 | Measured (Fig 6)                      |
| Dampings              | Known (full eigenvalue picture)       | Not reported                          |
| Flutter               | Published boundary                    | Subcritical fold, onset ~13.19 m/s    |
| If a model disagrees  | Code bug (unambiguous)                | Code, parameters, *or* physics        |


## How this maps to the fidelity ladder

The contribution of the dissertation is a cost-accuracy frontier across aerodynamic fidelities (quasi-steady → Peter's finite-state → UVLM → CFD). Each rung must be **verified correct on Isogai before it is validated on Michigan**. Only then can the gap between a model's prediction and the experiment be attributed to *missing physics* (e.g. quasi-steady omitting the wake-memory lag) rather than to an undetected bug. Verify first, validate second — otherwise the entire frontier sits on unverified code.

## Supervisor framing (10 June)

Fintan endorsed exactly this split: use a numerical benchmark with full frequency-and-damping data to confirm the numerics are correct, while the experiment (Michigan) remains the key comparison at the end. So Michigan is the headline *experimental* reference, with a numerical benchmark retained specifically for verifying numerics — both, for different jobs.

## "Match" means compare, not fit

A common slip is to read the workflow as "once Isogai checks the structure, aim the aero model at Michigan until it matches." That is the one thing not to do. *Match* here means **compare against**, never **tune until it agrees**. The aero model is never adjusted to hit the experiment.

Two clarifications make this precise:

- **Isogai verifies the aero too, not just the structure.** Verification is not "structure once, then the aero meets Michigan for the first time." Example 0 verifies the skeleton; the Stage 1 VGBF on Isogai verifies the *quasi-steady aero* against Isogai's known flutter behaviour; Peter's, UVLM, and CFD each get the same Isogai pass later. The aero is always verified on Isogai *before* it is compared on Michigan.
- **Only the structure is aimed at Michigan.** The Michigan calibration touches *structural / sectional* properties only — the inertia and stiffness ratios tuned to the measured 5.3 / 6.2 Hz and the rig's damping. The aero stays fixed.

The pipeline per fidelity rung is therefore:

1. **Verify on Isogai** — the structure (example 0) and that rung's aero (its VGBF against the known benchmark).
2. **Calibrate the Michigan structure** to the measured modes (5.3 / 6.2 Hz, damping). Aero fixed.
3. **Run the verified aero on the calibrated section** and *compare* the prediction to the experiment.
4. **The mismatch per rung is the data point.** No tuning of the aero, ever.

And the reassurance: quasi-steady *not* matching Michigan is the expected, correct outcome — it is the first frontier point, "here is how far off the cheapest model is." As the fidelity climbs to Peter's, UVLM, and CFD, the predicted lines should creep toward the experiment, and *that convergence is the result*. If each aero were instead tuned to match, every line would sit on the experimental curve and nothing about cost versus accuracy would have been demonstrated.