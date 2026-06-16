"""Null Object Pattern

Special case of no aerodynamic forces, no added aero states, aero rhs also 0 for stage 0."""

import numpy as np

class NoAero: #no need to implement protocol again since done in base.py
    """Null Object for NoAero aerodynamic models. Returns 0 states, (2,) vector for forces and empty array for (0,) aero rhs"""
    #for protocol, methods are stubs (...) but for implementing classes, instances return something
    @property
    def n_aero_states(self) -> int:
         """returns 0 for extra ODE terms"""
         return 0
    def forces(self, tau, y_struct, y_aero, U_star) -> np.ndarray:
         """returns zeroes of (2,) shape for aerodynamics forces"""
         return np.zeros(2)  #2 bc 2-DOF, pitch-plunge
    def aero_rhs(self, tau, y_struct, y_aero, U_star) -> np.ndarray:
         """returns zeroes of (n_aero_states,) shape for aero_rhs"""
         return np.zeros(0) #len(aero_rhs(...)) == n_aero_states
    def K_aero(self, U_star):
        return np.zeros((2, 2))
    def C_aero(self, U_star):
        return np.zeros((2, 2))