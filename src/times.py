"""
Timing equations.
"""
from numpy import cos

from constants import c
from lens import angle_eqn
from numeric import arr_roots

def time_curve(θ, cₘ, D_L, D_S, θ_S):
    return θ / (2*cₘ) / (1/(D_L * cos(θ_S)) - 1/D_S)

def time_stretch(Φ2D, cₘ):
    return 2*Φ2D / (cₘ * c**2)

def time_delay(Φ2D, diffΦ2D, cₘ, η, D_L, D_S, θ_, θ_S):
    eqn = angle_eqn(Φ2D, diffΦ2D, η, D_L, D_S, θ_, θ_S)

    roots = arr_roots(eqn)

    θ_images = θ_[1:-1][roots] - θ_S
    Φ2D_images = Φ2D[1:-1][roots] - θ_S

    T_c = time_curve(θ_images, cₘ, D_L, D_S, θ_S)
    T_s = time_stretch(Φ2D_images, cₘ)

    return T_c + T_s
