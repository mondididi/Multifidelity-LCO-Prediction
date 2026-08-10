# SU2 prescribed-pitch campaign — run notes (cmd)

Goal: the inviscid-ceiling figure (F8). Two forced-pitch cases at the rig's
reduced frequency k = 0.286: 2 deg (attached benchmark vs Peters) and 14 deg
(sweeps -4..+24 deg — the Euler loop must stay an attached ellipse straight
through the stall range: no hysteresis, no fold mechanism at any inviscid
fidelity). Confirmatory, not load-bearing: the paper stands without it.

## Files
Copy `naca0020_medium.su2` into this folder first (the cfgs inherit the
steady setup and expect the mesh beside them).

## Run (from CFD\SU2\Euler\unsteady, mflco env active)

    copy ..\naca0020_medium.su2 .
    SU2_CFD steady_a10.cfg
    copy restart_steady10.dat solution_flow.dat    & rem safety: some builds read this name
    SU2_CFD pitch_small_k0286.cfg
    SU2_CFD pitch_large_k0286.cfg
    set PYTHONPATH=..\..\..\..\src
    python su2_unsteady_energy.py

## What to expect
- steady_a10: minutes. Watch CL settle; at 10 deg inviscid CL ~ 1.1-1.2
  (2*pi*alpha-ish — NOT the stalled 0.76: that is the whole point).
- each unsteady case: 512 physical steps x 80 inner iters. Rough budget
  1-4 h/case on a laptop — launch overnight. Screen shows TIME_ITER climbing;
  inner RMS_DENSITY should drop ~2-3 orders each step.
- post-script: F8_euler_loops.png + su2_unsteady_energy.json here.
  Expected: small-case loop lying on the Peters loop; large-case a fat
  attached ellipse crossing the static-polar's stall plateau; cycle work
  W = closed-int Cm dalpha NEGATIVE (aero damping) for both.

## If something misbehaves
- "SOLUTION_FILENAME not found" on the unsteady cases: set RESTART_SOL= NO
  in the cfg (freestream start; first transient period absorbs it — that is
  why 4 periods are simulated and only the last 2 are used).
- Inner iterations not converging at the 14-deg case (transonic pockets:
  peak incidence 24 deg at M 0.3 can locally accelerate hard): halve
  TIME_STEP to 4.2011e-4 and double TIME_ITER to 1024. If it still fights,
  drop MACH_NUMBER to 0.2 in BOTH steady_a10 and pitch_large (k is what
  matters; re-derive PITCHING_OMEGA = 0.286*0.2*340.26/0.5 = 38.94 and
  TIME_STEP = 2*pi/38.94/128 = 1.2605e-3).
- Hard stop remains Tue 11 EOD: if the large case will not run, ship the
  small case alone (ceiling vs Peters) or steady-only + positioned, both
  pre-agreed as acceptable.
