#!/usr/bin/env python3
"""
Main interface.
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.constants import degree

from argparse import ArgumentParser

import lens
import potential
from util import mirror


def get_args():
    parser = ArgumentParser(prog="GravLens")

    io = parser.add_argument_group("IO")
    val = parser.add_argument_group("Values")

    io.add_argument("-o", "--output",
        required=True,
        help="Output directory.")

    val.add_argument("--radius", type=float,
        default=2.2,
        help="Radius of galaxy cluster in Mpc (default 2.2)")
    val.add_argument("--dist-lens", type=float,
        default=54.0,
        help="Distance to galaxy cluster in Mpc (default 54.0)")
    val.add_argument("--dist-source", type=float,
        default=10.0,
        help="Distance to source in units of distance to lens (default 10.0)")
    val.add_argument("--angle-source", type=float,
        default=10.0,
        help="Angle between source and lens in units of inverse angular "
             "radius of lens (default 10.0)")
    val.add_argument("--mass", type=float,
        default=1.2e15,
        help="Mass of galaxy cluster in Msun (default 1.2e15)")
    val.add_argument("--profile-index", type=float,
        default=1.0,
        help="Index on mass profile, i.e. 'γ' (default 1.0)")
    val.add_argument("--scaling-radius", type=float,
        default=1e-10,
        help="Scaling radius on mass profile, i.e. 'a' (default 1e-10)")

    args = parser.parse_args()

    args.dist_source *= args.dist_lens

    return args


def main():
    args = get_args()

    R = args.radius
    D_L = args.dist_lens
    D_S = args.dist_source
    M = args.mass

    γ = args.profile_index
    a = args.scaling_radius

    θmax = np.arctan2(R, D_L)
    θ = np.linspace(0, θmax, 1000)
    dθ = θ[1] - θ[0]

    θ_S = θmax / args.angle_source

    θ2 = mirror(θ, negate=True)
    θ_2 = θ2 + θ_S

    Φ2D = potential.flatten(potential.dehnen3D,
                            R, D_L, θ,
                            γ=γ, M=M, a=a)
    Φ2D2 = mirror(Φ2D)

    diffΦ2D2 = np.diff(Φ2D2) / dθ

    ωₚ = None
    ωₗ = None
    η = 1

    eqn = lens.angle_eqn(Φ2D2, diffΦ2D2, η, D_L, D_S, θ_2, θ_S)
    angles = lens.incident_angles(Φ2D2, diffΦ2D2, η, D_L, D_S, θ_2, θ_S)

    fig, ax1 = plt.subplots()

    ax1.plot(θ_2, Φ2D2, "k-")


    ax2 = ax1.twinx()

    ax2.plot(θ_2[:-1], diffΦ2D2, "r--")


    fig, ax = plt.subplots()

    ax.plot(θ_2[:-1], eqn)

    for θ_ in angles + θ_S:
        ax.axvline(θ_, color="red", linestyle="--")

    print(angles)
    plt.show()


    return


    plt.plot(θ/degree, Φ2D)
    plt.show()

    plt.plot(θ[1:]/degree, diffΦ2D)
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
