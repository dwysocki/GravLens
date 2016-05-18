"""
Density equations.
"""
import numpy as np
from scipy.integrate import quad

def average(density, R, D, θ, *args, **kwargs):
    ρ = np.empty_like(θ)
    d = D*np.tan(θ)
    s_max = np.sqrt(R**2 - d**2)

    def r(s, d):
        return np.sqrt(s**2 + d**2)

    for i in range(len(d)):
        def density_(s):
            return density(r(s, d[i]), *args, **kwargs)

        ρ[i] = 2*quad(density_, 0, s_max[i])[0]

    return ρ / (2 * s_max)



def dehnen3D(r, *, γ, M, a):
    const = (3 - γ) * M * a / (4*np.pi)
    rest  = r**-γ * (r+a)**(γ - 4)

    return const * rest
