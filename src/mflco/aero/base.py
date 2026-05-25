"""Define Aeromodel base class.

Solver needs something that provides forces/aero-state derivatives.
if a class contains n_aero_states,forces,aero_rhs, then it is classified as an aeromodel.
"Floor" of the aeromodel (minimum requirement) using protocol to ensure compatible data input/output for all aeromodels.
we define the Protocol because the solver requires a specific interface, and writing that interface down explicitly is better than leaving it implicit. """


from typing import Protocol
import numpy as np

class AeroModel(Protocol):
    """Protocol for aero model, which defines the required interface for any aero model to be used in the solver."""
    @property #add property decorator so when calling it doesn't need parentheses, it behaves like an attribute (bc its fixed descriptor)
    def n_aero_states(self) -> int:  #self is required as first parameter for method in a class, so it knows which obj to use
        """return how many extra ODE states the aero adds. Zero for NoAero and quasi-steady, a few for WJ and many for UVLM."""
        ...
    def forces(self,tau, y_struct, y_aero, U_star) -> np.ndarray:
        """returns 2 vector, Q_xi (plunge) and Q_alpha (pitch) , which are generalized aero forces for feeding into structural RHS."""
        ...
    def aero_rhs(self, tau, y_struct, y_aero, U_star) -> np.ndarray:
        """returns a vector of length n_aero_states for rhs""" 
        ...