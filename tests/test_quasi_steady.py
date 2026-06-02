"""Unit tests for the QuasiSteady aero model (circulatory-only, Prandtl-Glauert).

Tests forces() in isolation: each uses a hand-picked y_struct with a known
expected output, so a failure localises to the force model itself, not the
solver. Five tests, each framed by the bug class it catches:

1. Zero state -> zero force.
   Every term in `circ` is state-proportional, so a zero state must give
   exactly [0, 0]. Catches a stray constant/offset leaking into the force.

2. Steady thin-airfoil lift.
   Frozen pose (alpha only): circ -> U*^2 * alpha, so Q_xi -> -(2/mu)*U*^2*alpha
   at M=0. The generalized-force image of C_L -> 2*pi*alpha (the pi cancels
   against the 1/(pi*mu) prefactor). Catches a wrong prefactor: leftover pi,
   or 1/mu instead of 2/mu.

3. Prandtl-Glauert scaling.
   Same state at M=0 vs M=0.5 must differ purely by 1/sqrt(1-M^2) (~1.1547).
   Tested as a ratio, so it is independent of the prefactor (test 2).
   Catches comp_factor missing, misplaced, or applied with the wrong power.

4. Moment-lift invariant.
   Both forces share the same circulatory bracket, so Q_alpha = -(1/2+a)*Q_xi
   at every state. Checked at a fully nonzero state (all rates, M>0) so every
   term in `circ` and comp_factor is active. Catches a typo in one force's
   term that does not appear in the other.

5. Aero-state contract.
   QS declares zero memory states; aero_rhs must return a length-0 array so the
   solver's setup-time assert holds. The degenerate case of a contract that
   becomes load-bearing at Stage 2.

Reference (M=0 so comp_factor=1 unless noted):
    circ    = U*^2 * alpha + U* * xi_dot + U* * (1/2 - a) * alpha_dot
    Q_xi    = -(2/mu) * comp_factor * circ
    Q_alpha = +(2/mu) * (1/2 + a) * comp_factor * circ
"""

import numpy as np
import pytest
from mflco.model.params import TypicalSectionParameters
from mflco.aero.quasi_steady import QuasiSteady

p = TypicalSectionParameters() #default call, no beta addition, same case for all unit tests
a = p.a   #universal, geometric
mu = p.mu #universal, geometric

qs = QuasiSteady(p, M_inf = 0.0) #quasi steady at approximately M = 0, picks Mach t 0.0, only change in certain cases
    #this line creates a qs object that has 3 attributes, n_aero_states, forces, and aero_rhs
    #need to call out forces to get the Q_xi & Q_alpha


def test_zero_state_gives_zero_force():
    test_y_struct = [0.0, 0.0, 0.0, 0.0] #[xi, alpha, xi_dot, alpha_dot], all zeros bc zero_state
    test_U_star = 1.5

    Q = qs.forces(tau=0.0,y_struct=test_y_struct,y_aero=[],U_star=test_U_star) #but its unit test, so set y_struct and U_star as anything, tau and aero 0.0 & []

    assert Q[0] == 0.0
    assert Q[1] == 0.0  


def test_steady_lift_thin_airfoil():
    test_y_struct =[0.0, 0.1, 0.0, 0.0]
    test_U_star = 1.5

    Q = qs.forces(tau=0.0,y_struct=test_y_struct,y_aero=[],U_star=test_U_star)
    Q_xi_hand = -(2.0/mu) * (test_U_star**2) * (test_y_struct[1])
    Q_alpha_hand = (2.0/mu) * (0.5+a) * (test_U_star**2) *(test_y_struct[1])

    assert Q[0] == pytest.approx(Q_xi_hand)
    assert Q[1] == pytest.approx(Q_alpha_hand)


def test_prandtl_glauert_scaling():
    test_y_struct =[0.0, 0.1, 0.0, 0.0]
    test_U_star = 1.5

    qs_at_M5 = QuasiSteady(p, M_inf=0.5)
    Q_at_M5 = qs_at_M5.forces(tau=0.0,y_struct=test_y_struct,y_aero=[],U_star=test_U_star)
    Q_at_M0 = qs.forces(tau=0.0,y_struct=test_y_struct,y_aero=[],U_star=test_U_star)
    
    #check if differ by ~= 1.1547
    assert Q_at_M5[0] == pytest.approx(Q_at_M0[0]*(1/np.sqrt(1-0.5**2)))
    assert Q_at_M5[1] == pytest.approx(Q_at_M0[1]*(1/np.sqrt(1-0.5**2)))
    

def test_moment_lift_invariant():
    qs_case_4 = QuasiSteady(p,M_inf = 0.6) #M_inf > 0

    test_y_struct = [0.02, 0.1, -0.03, 0.05] #all non-zeroes
    test_U_star = 1.5

    Q_case_4 = qs_case_4.forces(tau=0.0,y_struct=test_y_struct,y_aero=[],U_star=test_U_star)
    Q_xi = Q_case_4[0]
    Q_alpha = Q_case_4[1]

    assert Q_alpha == pytest.approx(-(0.5+a)*Q_xi)


def test_aero_state_contract():
    test_y_struct =[0.0, 0.1, 0.0, 0.0]
    test_U_star = 1.5

    assert qs.n_aero_states == 0 #check if number of aero state eqs is 0
    assert qs.aero_rhs(tau=0.0,y_struct=test_y_struct,y_aero=[],U_star=test_U_star).shape == (qs.n_aero_states,) #check if return shape (0,)