import numpy as np
from mflco.model.michigan_params import calibrate_michigan, section_from_params
from mflco.model.analysis import undamped_natural_frequencies

cal = calibrate_michigan()
sec = section_from_params(cal)
nat, _ = undamped_natural_frequencies(sec)          # nondim omega/omega_alpha
print(np.sort(nat) * cal.omega_alpha / (2*np.pi))   # expect ~ [5.3, 6.2]