#refer to (1) Isogai (1979) "On the transonic-dip mechanism of flutter of a sweptback wing"
#         (2) García Pérez (2024) "Data-Driven Bifurcation Analysis of Experimental Aerolastic systems using preflutter measurements"
#         (3) Bristol (2019)
     
# stores isogai-style non-dim parameters and assembles symmetric lagrangian-form M, C, K

import numpy as np

class TypicalSectionParameters():   #class handle of structural param (based on Isogai A),
    def __init__(self, 
                 a = -2.0, 
                 x_alpha = 1.8,
                 r_alpha_sq = 3.48, 
                 omega_ratio = 1.0, 
                 mu = 60.0, 
                 beta = 0.0, 
                 zeta_h = 0.0, 
                 zeta_alpha = 0.0):
        self.a = a
        self.x_alpha = x_alpha
        self.r_alpha_sq = r_alpha_sq
        self.omega_ratio = omega_ratio
        self.mu = mu
        self.beta = beta
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
        # K = [omega_ratio*omega_ratio,         0]
        #     [0,   r_alpha_sq*(1 + beta*alpha*alpha)]
        K = np.array([[self.omega_ratio**2, 0], 
                      [0, self.r_alpha_sq * (1 + self.beta * alpha**2)]]) # != self.alpha, 0 by default
        return K
        