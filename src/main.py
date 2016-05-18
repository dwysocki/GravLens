#!/usr/bin/env python3
"""
Main interface.
"""
import numpy as np
from numpy import pi
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.constants import arcmin

from argparse import ArgumentParser
from os import path

from constants import c
import density
import lens
import potential
import times
from util import make_sure_path_exists, mirror, mirrorlr


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
        default=100.0,
        help="Angle between source and lens in units of inverse angular "
             "radius of lens (default 10.0)")
    val.add_argument("--mass", type=float,
        default=1.2e15,
        help="Mass of galaxy cluster in Msun (default 1.2e15)")
    val.add_argument("--profile-index", type=float,
        default=0.1,
        help="Index on mass profile, i.e. 'γ' (default 1.0)")
    val.add_argument("--scaling-radius", type=float,
        default=1e-2,
        help="Scaling radius on mass profile, i.e. 'a' (default 1e-10)")
    val.add_argument("--frequencies", type=float, nargs="+",
        default=[0.9, 1.0, 1.1],
        help="List of photon frequencies in GHz (default 0.9, 1.0, 1.1)")
    val.add_argument("--fraction-plasma", type=float,
        default=0.10,
        help="Fraction of total mass in plasma (default 0.10)")

    args = parser.parse_args()

    args.dist_source *= args.dist_lens
    args.frequencies = np.multiply(1e9, args.frequencies)

    make_sure_path_exists(args.output)

    return args


def main():
    args = get_args()

    output = args.output

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

    ρ2D = density.average(density.dehnen3D,
                          R, D_L, θ,
                          γ=γ, M=M, a=a)
    nₑ = lens.electron_number_density(ρ2D)

    ωₚ = lens.angfreq_plasma(ρ2D*args.fraction_plasma)
    ωₗ_arr = 2*pi * args.frequencies

    n_plasma = np.array([
        lens.n_plasma(ωₚ, ωₗ)
        for ωₗ in ωₗ_arr
    ])
    n_plasma2 = mirrorlr(n_plasma)


    ηGW = 1
    ηEM = (n_plasma2 + 1) / 2

    GWeqn = lens.angle_eqn(Φ2D2, diffΦ2D2, ηGW, D_L, D_S, θ_2, θ_S)
    GWangles = lens.image_angles(Φ2D2, diffΦ2D2, ηGW, D_L, D_S, θ_2, θ_S)
    GWtimes = times.time_delay(Φ2D2, diffΦ2D2, c, ηGW, D_L, D_S, θ_2, θ_S)

    EMeqn = lens.angle_eqn(Φ2D2, diffΦ2D2, ηEM[0,1:], D_L, D_S, θ_2, θ_S)
    EMangles = [
        lens.image_angles(Φ2D2, diffΦ2D2, η[1:], D_L, D_S, θ_2, θ_S)
        for η in ηEM
    ]
    EMtimes = [
        times.time_delay(Φ2D2, diffΦ2D2, c/n, η[1:], D_L, D_S, θ_2, θ_S)
        for n, η in zip(n_plasma2, ηEM)
    ]

    print("GW angle:", GWangles)
    print("GW time:", GWtimes)

    print("EM angle:", list(zip(args.frequencies/1e9, EMangles)))
    print("EM time:", list(zip(args.frequencies/1e9, EMtimes)))

    ##############
    ## Plotting ##
    ##############

    # Plot potential and its derivative.
    fig, ax_potential = plt.subplots()

    p_potential = ax_potential.plot(θ2/arcmin, Φ2D2, "k-",
                                    label=r"$\Phi_{\rm 2D}$")

    ax_diff_potential = ax_potential.twinx()
    p_diff_potential = ax_diff_potential.plot(θ2[:-1]/arcmin, diffΦ2D2*arcmin,
                                              "r--",
                                              label=r"$\Phi_{\rm 2D}'$")

    ax_potential.set_xlabel(r"$\theta$ [arcmin]")
    ax_potential.set_ylabel(r"$\Phi_{\rm 2D}$ [Mpc$^3$ s$^{-2}$]")
    ax_diff_potential.set_ylabel(
        r"$\Phi_{\rm 2D}'$ [Mpc$^3$ s$^{-2}$ arcmin$^{-1}$]")

    ps = p_potential + p_diff_potential
    labels = [p.get_label() for p in ps]

    ax_potential.legend(ps, labels, loc="lower right")

    fig.savefig(path.join(output, "potential.pdf"))

    plt.close(fig)
    del fig, ax_potential, ax_diff_potential, p_potential, p_diff_potential


    # Plot average mass density, and number density of electrons in plasma.
    fig, ax_mass_density = plt.subplots()

    p_mass_density = ax_mass_density.plot(θ2/arcmin, mirror(ρ2D), "k-",
                                          label=r"$\rho_{\rm 2D}$")

    ax_number_density = ax_mass_density.twinx()
    p_number_density = ax_number_density.plot(θ2/arcmin, mirror(nₑ), "g-",
                                              label=r"$n_{\rm e}$")

    ax_mass_density.semilogy()
    ax_number_density.semilogy()

    ax_mass_density.set_xlabel(r"$\theta$ [arcmin]")
    ax_mass_density.set_ylabel(r"$\rho_{\rm 2D}$ [M$_\odot$ Mpc$^{-3}$]")
    ax_number_density.set_ylabel(r"$n_{\rm e}$ [Mpc$^{-3}$]")

    ps = p_mass_density + p_number_density
    labels = [p.get_label() for p in ps]

    ax_mass_density.legend(ps, labels, loc="upper right")

    fig.savefig(path.join(output, "density.pdf"))
    plt.close(fig)
    del fig, ax_mass_density, ax_number_density


    # Plot image angles over the equation that finds them.
    fig, (axEM, axGW)  = plt.subplots(2, sharex=True)

    axGW.plot(θ2[:-1]/arcmin, GWeqn/arcmin)
    axEM.plot(θ2[:-1]/arcmin, EMeqn/arcmin)

    axGW.axhline(0, color="black", linestyle="--")
    axEM.axhline(0, color="black", linestyle="--")

    for θ in GWangles:
        axGW.axvline(θ/arcmin, color="red", linestyle="--")
    for θ in EMangles[0]:
        axEM.axvline(θ/arcmin, color="red", linestyle="--")

    axEM.set_xlabel(r"$\theta$ [arcmin]")
    axEM.set_ylabel(r"EM Equation [arcmin]")
    axGW.set_ylabel(r"GW Equation [arcmin]")

    fig.savefig(path.join(output, "image-angle-equation.pdf"))

    plt.close(fig)
    del fig, axEM, axGW


    # # Plot image angles for GW and EM, colored by arrival times
    # fig = plt.figure()
    # gs = mpl.gridspec.GridSpec(2, 1, height_ratios=[5, 1])
    # ax_EM = fig.add_subplot(gs[0])
    # ax_GW = fig.add_subplot(gs[1], sharex=ax_EM)

    # ax_GW.scatter(GWangles/arcmin, np.zeros_like(GWangles),
    #               color="black", marker="|")

    # ax_GW.set_yticks([])
    # ax_GW.set_xlabel(r"$\theta$ [arcmin]")
    # ax_EM.set_ylabel(r"$\omega_\ell$ [rad s$^{-1}$]")

    # fig.savefig(path.join(output, "image-angle-GW-EM.pdf"))

    # plt.close(fig)
    # del fig, ax_EM, ax_GW



if __name__ == "__main__":
    exit(main())
