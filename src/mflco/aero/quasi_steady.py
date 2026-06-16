"""Quasi-Steady aerodynamic model.

        Use Theodorsen's theory but with the function C(k) = 1, which means wake responds instantaneously
        Use Prandtl-Glauert Compressibility correct for subsonic flow

        Obtain [Q_xi,Q_alpha]' """

import numpy as np

class QuasiSteady:

    def __init__(self,params, M_inf): #define QuasiSteaady as a callable
        """Callable function QuasiStady(p,M_inf=) 
        where p = TypicalSectionParameters and
        M_inf is an input to be defined, fixed Mach
        
        returns an object with 3 attributes, n_aero_states,
        forces, and aero_rhs. For QS, only forces is non_stub.
        """
        #---
        """Constructor, attach M_inf to QS class
        
        callable function: qs = QuasiSteady(p, M_inf = 0.8)
        Defined bc the QuasiSteady needs to remember that it
        needs M_inf as an extra term while carrying the section
        params (reference)"""

        self.params = params
        self.M_inf = M_inf

    @property #only provide n_aero_states as a property, bc directly below @ prop
    def n_aero_states(self) -> int: #no additional aero euqations
        return 0 
    
    def forces(self, tau, y_struct, y_aero, U_star) -> np.ndarray:
        """A function to return Aero Forces (Q_xi, Q_alpha) in shape (2,)
        
        in QS, use compressibility factor (Prandtl-Glauert) and with
        Theodorsen's C(K) = 1 assumption. Also dropped apparent mass because
        (1) its a high-reduced-frequency effect, (2) its going to make EOM implicit. 
        The final form has the compressibiliy factor attached to the circulation term 
        (derived from alpha_effective)
        
        Tau is unused, and y_aero is unused. This is because y_aero is wake-memory states,
        and QS has none (because assumes wake responds instantaneously)"""
        [xi, alpha, xi_dot, alpha_dot] = y_struct

        mu = self.params.mu #always symbolic
        a = self.params.a   #always symbolic
        M_inf = self.M_inf 

        comp_factor = 1.0 / np.sqrt(1.0 - M_inf**2) #Prandtl-Glauert Compressibility factor

        #shared circulatory bracket, refer to 3.3 in aero_model. md
        circ = (U_star**2 * alpha) + (U_star * xi_dot) + (U_star * (0.5-a) * alpha_dot) #a is from self, based on geometry isogai
        #the circulatory term is the one derived from alpha_eff, so it must be multiplied with comp factor (3.3 aero_model)

        #aero forces
        Q_xi = (-2.0/mu)*circ*comp_factor #was -2/mu*circ, but added compressibility, 
        #-1/mu term is neglected bc contains ddot terms, which are apparent masses
        Q_alpha = (2.0/mu)*(0.5+a)*circ*comp_factor #refer to 3.3

        return np.asarray([Q_xi, Q_alpha]) #shape (2,), order plunge, pitch. 2 vector that matches eom.py [xi, alpha]

    def aero_rhs(self, tau, y_struct, y_aero, U_star) -> np.ndarray:    #wake memory terms: if no wake memory, nothing to integrate.
        return np.zeros(0)
    
    def K_aero(self, U_star):
        comp_factor = 1/np.sqrt(1 - self.M_inf**2)
        g = (2/self.params.mu)*comp_factor
        K_mat = np.asarray([[0, g*U_star**2],
                            [0, -g*(0.5 + self.params.a)*U_star**2]])
        return K_mat

    def C_aero(self, U_star):
        comp_factor = 1/np.sqrt(1 - self.M_inf**2)
        g = (2/self.params.mu)*comp_factor
        C_mat = np.asarray([[g*U_star,                       g*U_star*(0.5 - self.params.a)],
                            [-g*(0.5 + self.params.a)*U_star, -g*(0.5 + self.params.a)*U_star*(0.5 - self.params.a)]])
        return C_mat