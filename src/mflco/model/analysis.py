"""Eigen Analysis Function to Compute Natural Frequencies and Mode Shapes of the Coupled System.

This module utilizes the Q matrices from Cooper (Adding trivial identitiy expression) to exapnd 
2x2 matrices (from symmetric M, C, K) to 4x4 matrices to allow complex conjugates to occur.
The complex conjugates can happen from damping and decay rates.

The Q matrix incorporates M, C, K into a single 4×4 matrix.
"""

from scipy.linalg import eigh #hermitian eigenvalue solver for symmetric
import numpy as np
from .params import TypicalSectionParameters

def undamped_natural_frequencies(params: TypicalSectionParameters, alpha_eq: float = 0.0):
    '''Compute natural frequencies and mode shapes from mass and stiffness matrices.
    
    Only valid for undamped system (C=0) and no aero (wind-off), since it uses the 
    generalized eigenvalue problem K v = lambda M v.'''
    M = params.mass_matrix()
    K = params.stiffness_matrix(alpha_eq)
    eigvalues, eigvecs = eigh(K, M) #generalized form
    natural_frequencies = np.sqrt(eigvalues)
    return natural_frequencies, eigvecs

def system_matrix(params: TypicalSectionParameters, alpha_eq: float =0.0, U_star: float =0.0, aero=None): #default values
    """ Build 4x4 first-order system Matrix Q (A in Cooper) for eigenvalue analysis
    
    The 2nd-order EOM is Mq_dotdot + C q_dot + K q = Q_aero, where q = [xi, alpha]. is converted to
    the 1st-order form dy/dt = Q*y, where y = [q q_dot] = [xi, alpha, xi_dot, alpha_dot].
    
    Q = [ 0 I
          -M^-1 K  -M^-1 C ] 
          
    when aero is provided, M_aero, C_aero(U_star), K_aero(U_star, alpha_eq) are added to the structural matrices before assembly.

    Parameters:
     ----------
    params : TypicalSectionParameters
        Structural parameters.
    alpha_eq : float
        Pitch equilibrium angle (radians) about which to linearize K.
    U_star : float
        Nondimensional freestream speed. Zero by default (wind-off).
    aero : AeroModel or None
        Aero model providing linearized M_aero, C_aero, K_aero contributions.
        None for wind-off / structural-only analysis.
    
    Returns
    -------
    Q : np.ndarray, shape (4, 4)
        The first-order system matrix.
    """
    M = params.mass_matrix()
    C = params.damping_matrix()
    K = params.stiffness_matrix(alpha_eq) #evaluated at alpha_eq
    
    if aero is not None and  U_star > 0.0:
        #if aero is provided and U_star > 0, get aero contributions from protocol and add to structural matrices.
        raise NotImplementedError("Aero contributions to system matrix not implemented yet.")
        #M_aero, C_aero, K_aero
        #M = M + M_aero
        #C = C + C_aero
        #K = K + K_aero

    M_inv = np.linalg.inv(M) #for M^-1 in the Q matrix)

    Q = np.block([
        [np.zeros((2, 2)), np.eye(2)],     #identity of 2x2 bc [q, q_dot]
        [-M_inv @ K, -M_inv @ C] #@ is matrix multiplication
    ])

    return Q

def linearized_eigenvalues(params: TypicalSectionParameters, alpha_eq: float = 0.0, U_star: float = 0.0, aero=None):
    """Compute the 4 complex eigenvalues of the linearized system matrix Q, they represent the natural frequencies
    
    Builds Q via system_matrix, then computes np.linalg.eigvals(Q). The 4 eigenvalues
    organize as 2 complex-conjugate pairs, one pair per physical mode. Each pair encodes:
      - frequency      = |Im λ|
      - damping rate   = -Re λ  (negative real → decay, positive real → growth)

    Eigenvalues are sorted by |Im λ| descending, so high-frequency modes come first.

    Parameters
    ----------
    params, alpha_eq, U_star, aero : see system_matrix.

    Returns
    -------
    eigvalues : np.ndarray, shape (4,), dtype complex
        The 4 eigenvalues of Q, sorted by |Im λ| descending.

    **Only returns raw eigenvalues**
    """
    Q = system_matrix(params, alpha_eq, U_star, aero)
    eigvalues = np.linalg.eigvals(Q)    # np.linalg.eigvals handles the standard problem Q·v = λ·v (Q is not symmetric so we can't use eigh)

    #sort by |Im λ| descending, high f modes come first
    eigvalues = eigvalues[np.argsort(-np.abs(eigvalues.imag))]
    
    return eigvalues

def modal_analysis(params: TypicalSectionParameters, alpha_eq: float = 0.0, U_star: float = 0.0, aero=None):
    """Extract physical frequencies and damping ratios from the eigenvalues of the linearized system matrix Q.
    
    For each of the 2 physical modes, returns:
      - frequency ω = |Im λ|  (oscillation rate in units of ω_α)
      - damping ratio ζ = -Re λ / |λ|  (dimensionless decay rate)
        ζ < 0  → mode grows (unstable)
        ζ = 0  → mode oscillates without decay (neutral)
        ζ > 0  → mode decays (stable)
    
    Modes are returned sorted by frequency ascending (low-frequency mode first).
    
    Parameters
    ----------
    params, alpha_eq, U_star, aero : see system_matrix.
    
    Returns
    -------
    frequencies : np.ndarray, shape (2,)
        Natural frequencies of each mode (rad / nondimensional time).
    damping_ratios : np.ndarray, shape (2,)
        Damping ratio of each mode (dimensionless).
    """
    eigvalues = linearized_eigenvalues(params, alpha_eq, U_star, aero)

    #take one representative per complex-conjugate pair (+imag part)
    positive_eigvalues = eigvalues[eigvalues.imag > 0] #4 eigenvalues, 2 duplicates, +,- so only take 1 representative
    
    #edge-case: if purely reaal eigens, divergence occurs, static instability ***need to address when happen.
    # Expected: 2 oscillatory modes (one per DOF)
    # If less, we've encountered a non-oscillatory eigenvalue (divergence)
    if len(positive_eigvalues) != 2:
        raise ValueError(
            f"Expected 2 oscillatory modes but found {len(positive_eigvalues)}. "
            f"This indicates a non-oscillatory eigenvalue (likely divergence). "
            f"Full eigenvalues: {eigvalues}. Handle divergence case explicitly."
    )

    #extract physical quantitites
    frequencies = np.abs(positive_eigvalues.imag) #frequency = |Im λ|
    damping_ratios = -positive_eigvalues.real / np.abs(positive_eigvalues)

    #sort ascendeing
    order = np.argsort(frequencies)

    return  frequencies[order], damping_ratios[order]