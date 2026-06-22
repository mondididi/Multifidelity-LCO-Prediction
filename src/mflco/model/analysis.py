"""Eigen Analysis Function to Compute Natural Frequencies and Mode Shapes of the Coupled System.

This module utilizes the Q matrices from Cooper (Adding trivial identitiy expression) to exapnd 
2x2 matrices (from symmetric M, C, K) to 4x4 matrices to allow complex conjugates to occur.
The complex conjugates can happen from damping and decay rates.

The Q matrix incorporates M, C, K into a single 4×4 matrix.
"""

from scipy.linalg import eigh #hermitian eigenvalue solver for symmetric
from scipy.linalg import eig as generalized_eig #generalized A x = lambda E x (Peters descriptor)
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
        C = C + aero.C_aero(U_star)
        K = K + aero.K_aero(U_star)

    M_inv = np.linalg.inv(M) #for M^-1 in the Q matrix)

    Q = np.block([
        [np.zeros((2, 2)), np.eye(2)],     #identity of 2x2 bc [q, q_dot]
        [-M_inv @ K, -M_inv @ C] #@ is matrix multiplication
    ])

    return Q

def descriptor_matrices(params: TypicalSectionParameters, alpha_eq: float = 0.0,
                        U_star: float = 0.0, aero=None):
    """Build the (4+N)x(4+N) descriptor pair (E, A) for an aero model with inflow states.

    First-order generalized form  E . y' = A . y,  y = [xi, alpha, xi', alpha', Lambda]
    (N = aero.n_aero_states inflow states). Solve as the generalized eigenproblem
    A x = lambda E x (E is nonsingular here: det E = det(M_s - M_a) * det(A_bar) != 0).

        E = [[ I,            0,                  0      ],
             [ 0,            M_s - M_a,          0      ],
             [ 0,            aero_forcing,       A_bar  ]]

        A = [[ 0,            I,                  0            ],
             [ -(K_s+K_aero),-(C_s+C_aero+C_a),  K_aero_lambda],
             [ 0,            aero_forcing_vel,   -U* . I      ]]

    Reuses the QS circulatory K_aero/C_aero unchanged; the Peters additions are
    the apparent mass/damping (M_a, C_a), the inflow->structure coupling
    (K_aero_lambda), and the structure->inflow forcing (aero_forcing on q'',
    aero_forcing_vel on q'). Unlike system_matrix this is NOT gated on U_star>0:
    the aero blocks carry their own U* factors and vanish where they should, so
    the U*=0 oracle (which needs M_a and aero_forcing present) is handled here.

    Parameters
    ----------
    params, alpha_eq, U_star : see system_matrix.
    aero : AeroModel with inflow states (n_aero_states > 0); must provide
        K_aero, C_aero, M_a, C_a, K_aero_lambda, aero_forcing, aero_forcing_vel,
        and the A_bar attribute.

    Returns
    -------
    E, A : np.ndarray, each (4+N, 4+N).
    """
    N = aero.n_aero_states

    M_s = params.mass_matrix()
    C_s = params.damping_matrix()
    K_s = params.stiffness_matrix(alpha_eq)

    M_eff = M_s - aero.M_a()
    C_eff = C_s + aero.C_aero(U_star) + aero.C_a(U_star)
    K_eff = K_s + aero.K_aero(U_star)

    K_al = aero.K_aero_lambda(U_star)        # (2, N)
    A_bar = aero.A_bar                        # (N, N)
    F_acc = aero.aero_forcing()               # (N, 2)  c_bar (x) S       on q''
    F_vel = aero.aero_forcing_vel(U_star)     # (N, 2)  c_bar (x) [0,-U*] on q'

    I2 = np.eye(2)
    Z2 = np.zeros((2, 2))
    Z2N = np.zeros((2, N))
    ZN2 = np.zeros((N, 2))

    E = np.block([
        [I2,  Z2,    Z2N],
        [Z2,  M_eff, Z2N],
        [ZN2, F_acc, A_bar],
    ])
    A = np.block([
        [Z2,     I2,     Z2N],
        [-K_eff, -C_eff, K_al],
        [ZN2,    F_vel,  -U_star * np.eye(N)],
    ])
    return E, A

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
    if aero is not None and getattr(aero, "n_aero_states", 0) > 0:
        # Peters (or any inflow model): generalized descriptor eigenproblem
        E, A = descriptor_matrices(params, alpha_eq, U_star, aero)
        eigvalues = generalized_eig(A, E, right=False)   # solves A x = lambda E x
    else:
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
    # one representative per complex-conjugate pair. Threshold (not >0) so a real
    # lag eigenvalue carrying ~1e-15 numerical imag is not miscounted as oscillatory.
    oscillatory = eigvalues[eigvalues.imag > 1e-8]

    # divergence guard: a real eigenvalue with Re>0 is static divergence, which
    # this oscillatory reduction would otherwise hide. (Flutter is an oscillatory
    # mode crossing Re=0, kept below.) Won't fire on Isogai/Michigan.
    real_eigs = eigvalues[np.abs(eigvalues.imag) <= 1e-8]
    if np.any(real_eigs.real > 1e-6):
        raise ValueError(
            f"Divergence: real eigenvalue with positive real part. {eigvalues}"
        )

    # The 2 STRUCTURAL modes are the 2 *least-damped* oscillatory modes. With
    # Peters the N finite-state lag poles need not be real -- for N=3, eig(A_bar^-1)
    # has a complex pair, giving a heavily-damped oscillatory lag mode (zeta~0.88,
    # U*-independent). Selecting by smallest damping ratio isolates the lightly
    # damped structure (zeta~0.01-0.02, -> 0 at flutter) from the damped wake lag,
    # robustly across N (lag is a fast, damped process: always more damped than
    # the structural modes). For QS (no lag) there are exactly 2 oscillatory and
    # this is a no-op.
    if len(oscillatory) < 2:
        raise ValueError(
            f"Fewer than 2 oscillatory modes found ({len(oscillatory)}). {eigvalues}"
        )
    zeta_all = -oscillatory.real / np.abs(oscillatory)
    structural = oscillatory[np.argsort(zeta_all)[:2]]

    #extract physical quantitites
    frequencies = np.abs(structural.imag) #frequency = |Im λ|
    damping_ratios = -structural.real / np.abs(structural)

    #sort ascendeing
    order = np.argsort(frequencies)

    return  frequencies[order], damping_ratios[order]