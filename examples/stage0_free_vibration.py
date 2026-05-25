"""Example of Stage 0: Free Vibration"""

from mflco.model.params import TypicalSectionParameters
from mflco.model.solver import integrate
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np

p = TypicalSectionParameters() #default as Isogai, no overrides
y0 = [0.0, 0.05, 0.0, 0.0]  #small pitch pertubation
tau_span = (0.0, 30.0)      #30 dimensionl time units 
t_eval = np.linspace(*tau_span,3001)    #3001 timesteps
sol = integrate(p, y0, tau_span=tau_span, t_eval=t_eval)   #let default tolerance, method, and zero aero (default none), output sol.t, sol.y

#unpack
[xi, alpha, xi_dot, alpha_dot] = sol.y
tau = sol.t

#plot xi v tau, alpha v tau, on stacked subplots
fig, axs = plt.subplots(2, 1, figsize=(10,6))

for ax in axs:
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.xaxis.set_minor_locator(MultipleLocator(0.5))
    ax.grid(which='major', alpha=0.6)
    ax.grid(which='minor', alpha=0.2)
    ax.axhline(0, color='k', linewidth=0.75, alpha=0.5)
    ax.set_xlim(*tau_span)

axs[0].plot(tau,xi)
axs[0].set(xlabel = 'Tau (-)', ylabel='Xi (-)', title='Non-Dimensional Time vs Plunge ξ(τ)')

axs[1].plot(tau,alpha)
axs[1].set(xlabel = 'Tau (-)', ylabel='Alpha (rad)', title='Non-Dimensional Time vs Pitch α(τ)')

fig.suptitle('Stage 0: Free Vibration, Isogai Case A,  wind-off')
fig.tight_layout()
plt.show()