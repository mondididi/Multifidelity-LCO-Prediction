# Structural Model Reference

A standalone reference for the 2-DOF pitch–plunge structural model implemented
in `mflco.model`. Covers the equations of motion, the structural matrices, the
nondimensionalisation, and why each matrix has the form it does.

## References

The dimensional EOM and the canonical Isogai parameters are reproduced verbatim
in several places — easiest in the project repo to point at:

- Isogai (1979) — original Case A specification. AIAA J., 17(7), 793–795.
- Lee, Price & Wong (1999) — *Prog. Aero. Sci.* 35, 205–334. Authoritative
  treatment of the cubic-spring LCO formulation; gives post-Hopf bifurcation
  data for various μ, x_α, r_α², β_α.
- Yuan, Sandhu & Poirel (2021) — *J. Aerospace Eng.* 34(2): 04020117.
  Reproduces Isogai's Case A parameters cleanly and includes the wind-off
  coupled frequencies (0.713 ω_α, 5.34 ω_α).
- García Pérez et al. (2024) — *AIAA J.* 62(5), 1906–1915. Experimental 2-DOF
  airfoil with cubic spring; potential validation target alongside or instead
  of Bristol/Tartaruga (2019).