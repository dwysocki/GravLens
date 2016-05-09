"""
Gravitational potential equations.
"""
import numpy as np


def flatten(potential, R, D, θ, *args, **kwargs):
    def r(s, θ):
        return np.sqrt(s**2 + (D*np.tan(θ))**2)

    s = np.linspace(-R, R, 1000)
    return potential(r(s, θ), *args, **kwargs)



def dehnen3D(r, *, γ, M, a):
    power = 2 - γ
    const = -G*M / (a * power)
    depend = 1 - (r / (r+a))**power

    return const * depend




def zero2D(θ):
    return 0

def zero3D(x):
    return 0
