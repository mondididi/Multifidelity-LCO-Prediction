## 5. Energy of the system

For the symmetric Lagrangian form:


**Kinetic energy:**

```
KE = ½ q̇ᵀ M q̇
   = ½ (ξ'² + 2 x_α ξ' α' + r_α² α'²)
```


**Potential energy** (linear plunge spring + nonlinear pitch spring):

```
PE = ½ (ω_h/ω_α)² ξ²              (plunge)
   + ½ r_α² α² + ¼ r_α² β α⁴      (pitch, including cubic term)
```

The cubic-spring potential comes from integrating the pitch restoring force:

```
∫ k_α,0 (1 + β α²) α dα = k_α,0 (½ α² + ¼ β α⁴)
```

and the leading `r_α²` factor again comes from the row rescaling. Note: the
nondimensionalisation absorbs `k_α,0` into the unit choice (K_α/(I_α ω_α²) = 1).


**Total energy** E = KE + PE is conserved for undamped (ζ = 0) free vibration
(Q_aero = 0). The test in `test_eom_nonlinear.py` verifies this conservation to
relative drift < 10⁻⁶ over 20 nondimensional time units.
