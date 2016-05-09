"""
Lensing equations.
"""

import numpy as np
from numpy import linalg as la
from scipy.constants import c


def phase_delay(θ, ω, dₒₗ, dₗₛ, Φ2D):
    # Contribution due to geometry.
    geom = (dₒₗ + dₗₛ) * dₒₗ / (2 * dₗₛ * c) * la.norm(θ)**2
    # Contribution due to potential.
    potential = 2 * Φ2D(θ) * c**-3

    return ω * (geom - potential)
