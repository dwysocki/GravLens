import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.constants import degree, parsec, kilo, mega
from astropy.constants import M_sun


import lens
import potential


def main():
    R = 2.2
    D = 54
    M = 1.2e15
    θmax = np.arctan2(R, D)
    θ = np.linspace(0, θmax, 1000)

    Φ2D = potential.flatten(potential.dehnen3D,
                            R, D, θ,
                            γ=1, M=M, a=1)

    plt.plot(θ/degree, Φ2D)
    plt.show()

    return

    r = np.linspace(0, R, 1000)
    a = 1

    Φ1   = potential.dehnen3D(r, γ=1, M=M, a=a)
    Φ1_5 = potential.dehnen3D(r, γ=1.5, M=M, a=a)
    Φ2_5 = potential.dehnen3D(r, γ=2.5, M=M, a=a)

    plt.plot(r, -Φ1)
    plt.plot(r, -Φ1_5)
    plt.plot(r, -Φ2_5)
    plt.loglog()
    plt.show()

    return

    dₒₗ = dₗₛ = ω = 1
    Φ2D = potential.zero2D

    θx = θy = np.arange(-np.pi/4, np.pi/4, 0.01)
    θX, θY = np.meshgrid(θx, θy)

    Δφ = np.array([
        [lens.phase_delay(np.array([θx_, θy_]), ω, dₒₗ, dₗₛ, Φ2D)
         for θy_ in θy]
        for θx_ in θx
    ])

    plt.contour(θX, θY, Δφ)
    plt.show()

if __name__ == "__main__":
    exit(main())
