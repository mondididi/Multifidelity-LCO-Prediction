#eigenanalysis function to verify the stability of the system

from scipy.linalg import eigh #hermitian eigenvalue solver for symmetric
import numpy as np
from .params import TypicalSectionParameters

def coupled_natural_frequencies(params: TypicalSectionParameters):
    '''Compute natural frequencies and mode shapes from mass and stiffness matrices.'''
    M = params.mass_matrix()
    K = params.stiffness_matrix()
    eigvals, eigvecs = eigh(K, M) #generalized form
    natural_frequencies = np.sqrt(eigvals)
    return natural_frequencies, eigvecs