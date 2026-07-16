"""Parameter container for the 2-DOF pitch-plunge typical section.

Defines TypicalSectionParameters: holds Isogai-style nondimensional parameters
and assembles the symmetric Lagrangian-form M, C, K matrices.
Defaults reproduce Isogai (1979) Case A.

See docs/structural_model for derivation, sign conventions, and references.
"""

import numpy as np

class TypicalSectionParameters():   #class handle of structural param (based on Isogai A),
    def __init__(self, 
                 a = -2.0, 
                 x_alpha = 1.8,     #semi-chords
                 r_alpha_sq = 3.48, 
                 omega_ratio = 1.0, 
                 mu = 60.0, 
                 beta = 0.0, 
                 delta = 0.0,   # spring neutral offset from the operating point [rad].
                                # 0 = Woolston symmetric special case (pure cubic) = Stages 0-2.
                                # Nonzero regenerates the quadratic term of Lee (1999) Eq. 31.
                 zeta_h = 0.0, 
                 zeta_alpha = 0.0):
        self.a = a
        self.x_alpha = x_alpha
        self.r_alpha_sq = r_alpha_sq
        self.omega_ratio = omega_ratio
        self.mu = mu
        self.beta = beta
        self.delta = delta
        self.zeta_h = zeta_h
        self.zeta_alpha = zeta_alpha


    def mass_matrix(self): 
        # return mass matrix, 2x2
        # M = [1,         x_alpha   ]
        #     [x_alpha,   r_alpha_sq]
        M = np.array([[1, self.x_alpha], 
                      [self.x_alpha, self.r_alpha_sq]])
        return M


    def damping_matrix(self): 
        # return damping matrix, 2x2
        # C = [2*zeta_h*omega_ratio,     0                      ]
        #     [0,                        2*zeta_alpha*r_alpha_sq]
        C = np.array([[2 * self.zeta_h * self.omega_ratio, 0], 
                      [0, 2 * self.zeta_alpha * self.r_alpha_sq]])
        return C


    def stiffness_matrix(self, alpha=0.0): #if alpha is zero, the stiffness matrix will be linear
        # return stiffness matrix, 2x2
        # Lee et al. (1999) Prog. Aero. Sci. 35 Eq. 31: M(a) = b0 + b1*a + b2*a^2
        # + b3*a^3. b0 is balanced at equilibrium and drops out. Expanding a cubic
        # spring whose neutral point sits delta from the operating point gives
        # b2 = 3*delta*b3, so b2 is NOT free -- delta is bounded by the rig trim.
        #   b1 = r_alpha_sq  -- held FIXED, so flutter does NOT move with delta
        #   b2 = r_alpha_sq * 3*beta*delta
        #   b3 = r_alpha_sq * beta
        # delta = 0 recovers r_alpha_sq*(1 + beta*alpha^2) BIT-IDENTICALLY.
        # K = [omega_ratio*omega_ratio,   0]
        #     [0,   r_alpha_sq*(1 + 3*beta*delta*alpha + beta*alpha^2)]
        K = np.array([[self.omega_ratio**2, 0], 
                      [0, self.r_alpha_sq * (1
                                             + 3 * self.beta * self.delta * alpha
                                             + self.beta * alpha**2)]]) # != self.alpha, 0 by default
        return K