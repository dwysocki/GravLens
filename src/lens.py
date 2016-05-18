"""
Lensing equations.
"""

import numpy as np
from numpy import linalg as la

from constants import c, e, mₑ, mₚ, ε_0, μ_0
from numeric import arr_roots


def angfreq_plasma(nₑ):
    """
    Calculates the angular frequency of a plasma, `ωₚ`, with the given number
    density of electrons, `nₑ`.
    """
    return np.sqrt((nₑ * e**2) / (ε_0 * mₑ))


def permittivity_plasma(ωₚ, ωₗ):
    """
    Calculates the relative permittivity of the plasma, `ε`, given the
    frequency of the plasma, `ωₚ`, and of the light propogating through it,
    `ωₗ`.
    """
    return 1 - (ωₚ/ωₗ)**2


def speed_of_light_medium(μ, ε):
    """
    Calculates the speed of light in a medium, with relative permeability, `μ`,
    and relative permittivity, `ε`.
    """
    return np.sqrt(μ * ε)


def n_grav(Φ):
    """
    Calculates the gravitational index of refraction, `n`, for a given
    potential, `Φ`.
    """
    return 1 - 2*Φ / c**2


def n_plasma(ωₚ, ωₗ):
    """
    Calculates the plasma index of refraction, `n`, for a given plasma
    frequency, `ωₚ`, and light frequency, `ωₗ`.
    """
    μ = μ_0
    ε = permittivity_plasma(ωₚ, ωₗ)
    cₚ = speed_of_light_medium(μ, ε)

    return c / cₚ


def angle_eqn(Φ2D, diffΦ2D, η, D_L, D_S, θ_, θ_S):
    return (
        (2*η / c**2)
      * (1/(D_L * np.cos(θ_S)) - 1/D_S)
      * diffΦ2D
      - θ_[:-1]
    )


def incident_angles(Φ2D, diffΦ2D, η, D_L, D_S, θ_, θ_S):
    eqn = angle_eqn(Φ2D, diffΦ2D, η, D_L, D_S, θ_, θ_S)

    roots = arr_roots(eqn)

    return θ_[1:-1][roots] - θ_S


def phase_delay(θ, ω, dₒₗ, dₗₛ, Φ2D):
    # Contribution due to geometry.
    geom = (dₒₗ + dₗₛ) * dₒₗ / (2 * dₗₛ * c) * la.norm(θ)**2
    # Contribution due to potential.
    potential = 2 * Φ2D(θ) * c**-3

    return ω * (geom - potential)
