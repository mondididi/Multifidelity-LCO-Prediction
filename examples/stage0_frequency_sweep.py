"""Stage 0 validation: frequency sweep vs alpha_eq with cubic spring (beta sweep).

Each row is one beta value. Columns are the two modes (plunge-dominated, pitch-dominated).
Validates that the cubic stiffening beta·alpha² propagates through the eigenvalue analysis,
and that the magnitude of the effect scales with beta.
"""

from mflco.model.params import TypicalSectionParameters
from mflco.model.analysis import modal_analysis
import matplotlib.pyplot as plt
import numpy as np

betas = [3.0, 6.0]
alpha_eq_values = np.linspace(0.0, 0.3, 31) #31 pitch increments from 0.0 to 0.3 rad
# alpha_eq is the pitch angle around which we linearize the nonlinear stiffness K(alpha).
# Eigenvalue analysis is a linear tool, so we evaluate K at alpha_eq (treating it as fixed)
# and study small perturbations around that point.
# Sweeping alpha_eq tests whether the cubic stiffening beta·alpha² is correctly captured.

fig, ax = plt.subplots(1, 2, figsize=(12, 5), sharex=True)  #1 row, 2 cols, one per mode

for count, beta in enumerate(betas):
    p = TypicalSectionParameters(beta=beta)  #non-linear

    freqs_low = np.zeros_like(alpha_eq_values)
    freqs_high = np.zeros_like(alpha_eq_values) #same shape as alpha_eq

    for i, alpha_eq in enumerate(alpha_eq_values): #loop throughout all alpha eq values
        freqs, _ = modal_analysis(p, alpha_eq=alpha_eq) #output returns frequency and damping ratio, ignore the latter. 
        #freqeuency is sorted, so low f mode is [0], plunge dominated, high f mode is [1], pitch dominated
        freqs_low[i] = freqs[0] #plunge-dom
        freqs_high[i] = freqs[1] #pitch-dom

    ax[0].plot(alpha_eq_values, freqs_low, 'o-', label=f'beta = {beta}')
    ax[0].set_title('Plunge-dominated mode (omega1)')

    ax[1].plot(alpha_eq_values, freqs_high, 's-', label=f'beta = {beta}')
    ax[1].set_title('Pitch-dominated mode (omega2)')

for axis in ax.flat:    #.flat to iterate over the 2D axes array
    axis.set_xlabel('Equilibrium pitch angle alpha_eq (rad)')
    axis.set_ylabel('Modal frequency / omega_alpha (-)')
    axis.grid(True, alpha=0.3)
    axis.set_xlim(alpha_eq_values[0],alpha_eq_values[-1])   #-1 = last element
    axis.legend()

fig.suptitle('Stage 0 validation: frequency dependence on alpha_eq (beta sweep)')
fig.tight_layout()
plt.show()