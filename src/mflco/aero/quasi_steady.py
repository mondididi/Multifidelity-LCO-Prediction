"""Quasi-Steady aerodynamic model.

Use Theodorsen's function but with C(k) = 1, which means wake is instanteous achieve steady
Use Prandtl-Glaubert Compressibility correct for subsonic flow

Obtain [Q_xi,Q_alpha]' """

import numpy as np

class QuasiSteady:
    @property
    def n_aero_states(self) -> int:
        return
    
    def forces(self, tau, y_struct, y_aero, U_star) -> np.ndarray:
        return

    def aero_rhs(self, tau, y_struct, y_aero, U_star) -> np.ndarray:
        return
    