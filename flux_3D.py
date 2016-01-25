#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# flux_3D.py
#
# Copyright (C) 2016 Ivo Alxneit
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#



import numpy as np
import sys


def bin(x, y, z, P_per_area, nZ, nTheta, dZ, dTheta):
    binned = np.zeros((nZ, nTheta + 1))
    for (X, Y, Z) in zip(x, y, z):
        theta = np.pi + np.arctan2(Y, X)   # theta 0 ..2pi
        i = int(np.trunc(Z / dZ))
        j = int(np.trunc(theta / dTheta))
        # precautions for round off errors
        if i >= nZ:
            i = i - 1
        if j > nTheta:
            j = j - 1
        binned[i, j] += P_per_area[i]

    return binned


def bin_cone(n_z, H, n_theta, r):
    """
    Read output of 'get_flux' from stdin.
    Perform binning of the data on the grid 'n_z' by 'n_theta+1'.

        r: (radius_at_base,radius_at_top)
        H: (0, height_of_cone)

    return the binned data (flux density)
    """
    (x, y, z) = np.loadtxt(sys.stdin, usecols=(0,1,2), unpack=True)

    Dz = H[1] / n_z    # intervall in z direction
    Dtheta = 2 * np.pi / (n_theta+1)
    R = r[0] - r[1]
    dl = np.sqrt(H[1]**2 + R**2)   # length of line on surface of cone
    R_xy = np.linspace(r[0] + Dz * R / (2 * H[1]),
                       r[0] + Dz * R / (2 * H[1]) - n_z * Dz * R / H[1],
                       n_z)

    area = R_xy * Dtheta * dl
    binned = bin(x, y, z, 1/area, n_z, n_theta, Dz, Dtheta)

    return binned


def bin_cylinder(n_z, H, n_theta):
    """
    Read output of 'get_flux' from stdin.
    Perform binning of the data on the grid 'n_z' by 'n_theta+1'.

        H: (0, height_of_cylinder)

    return the binned data (flux density) and
    the scalar 'R' (radius of cylinder)
    """
    (x, y, z) = np.loadtxt(sys.stdin, usecols=(0,1,2), unpack=True)

    Dz = H[1] / n_z
    Dtheta = 2 * np.pi / (n_theta + 1)
    R = np.sqrt(x[0]**2 + y[0]**2)

    area = np.asarray([ R * Dtheta * Dz for i in range(n_z) ])
    binned = bin(x, y, z, 1/area, n_z, n_theta, Dz, Dtheta)

    return binned, R


def bin_sphere(n_z, H, n_theta):
    """
    Read output of 'get_flux' from stdin.
    Perform binning of the data on the grid 'n_z' by 'n_theta+1'.

        H: (min_height, max_height)

    return the binned data (flux density)
    """
    (x, y, z) = np.loadtxt(sys.stdin, usecols=(0,1,2), unpack=True)

    Dz = (H[1]-H[0]) / n_z
    Dtheta = 2 * np.pi / (n_theta + 1)
    R = np.sqrt(x[0]**2 + y[0]**2 + z[0]**2)
    h = np.linspace(H[0], H[0] + (n_z) * Dz, n_z + 1) / R
    phi1 = np.arccos(h[:n_z])
    Dphi = np.arccos(h[1:]) - phi1
    R_xy = R * np.sin(phi1)

    area = R**2 * Dphi * Dtheta
    binned = bin(x, y, z, 1/area, n_z, n_theta, Dz, Dtheta)

    return binned, R_xy


def R2xy(R, theta):
    """
    converts 'R', 'theta' to 'x','y'
    'R' can be scalar or vector of the same size as axis0 of 'theta'
    """
    x = R * np.sin(theta)
    y = R * np.cos(theta)

    return x, y


def grid_cone(nZ, H, nTheta, r):
    z = np.linspace(0, H[1], nZ)
    Dz = H[1] / nZ
    R = r[0] - r[1]
    theta = np.linspace(0, 2 * np.pi, nTheta + 1)

    radius = np.linspace(r[0] + Dz * R / (2 * H[1]),
                         r[0] + Dz * R / (2 * H[1]) - nZ * Dz * R / H[1],
                         nZ)
    Ri = np.transpose(np.asarray([radius for i in range(nTheta + 1)]))

    thetai, zi = np.meshgrid(theta, z)
    xi, yi = R2xy(Ri, thetai)

    return xi, yi, zi


def grid_cylinder(R, nZ, H, nTheta):
    z = np.linspace(0, H[1], nZ)
    theta = np.linspace(0, 2 * np.pi, nTheta + 1)

    (thetai, zi) = np.meshgrid(theta, z)
    (xi, yi) = R2xy(R, thetai)

    return xi, yi, zi


def grid_sphere(R, nZ, H, nTheta):
    z = np.linspace(H[0], H[1], nZ)
    theta = np.linspace(0, 2 * np.pi, nTheta + 1)
    Ri = np.transpose(np.asarray([R for i in range(nTheta + 1)]))

    (thetai, zi) = np.meshgrid(theta, z)
    (xi, yi) = R2xy(Ri, thetai)

    return xi, yi, zi


def plot_4D(flux, Limits, x, y, z):
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm, colors
    import matplotlib.pyplot as plt

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    if Limits is not None:
        norm = colors.Normalize(vmin=Limits[0], vmax=Limits[1])
    else:
        norm = colors.Normalize()

    s_m = cm.ScalarMappable(cmap=cm.jet, norm=norm)
    s_m.set_array(flux)

    flux_colors = cm.jet(norm(flux))
    surf = ax.plot_surface(x, y, z, rstride=1, cstride=1,
                           facecolors=flux_colors,
                           linewidth=0, antialiased=False,
                           shade=False)
    cbar = fig.colorbar(s_m)
    cbar.set_label('Flux')

    plt.show()

if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser(description='Plot data returned from "get_flux".')
    parser.add_argument('-z', '--nZ', metavar='N', type=int, default=10,
                        help='Number of intervals in z-direction')
    parser.add_argument('-a', '--nTheta', metavar='N', type=int, default=36,
                        help='Number of intervals in theta-direction')
    parser.add_argument('-Z', '--Zlimits', metavar='Z', nargs=2, type=float,
                        required=True,
                        help='Maximum z value of object')
    parser.add_argument('-F', '--Fluxlimits', metavar='F', type=float,
                        nargs=2,
                        help='Minimum and maximum of flux (for colorbar)')
    parser.add_argument('-t', '--type',
                        choices=['cone', 'cylinder', 'sphere'], required=True,
                        help='Type of body')
    parser.add_argument('-r', '--rlimits', metavar='R', nargs=2, type=float,
                        help='Minimum and maximum value of cone radius')
    parser.add_argument('-P', '--Pfactor', metavar='P', type=float,
                        default=1.0,
                        help='Scaling factor for flux "P_factor"')

    args = parser.parse_args()

    if args.type == 'cylinder':
        (flux, R) = bin_cylinder(args.nZ, args.Zlimits, args.nTheta)
        (x, y, z) = grid_cylinder(R, args.nZ, args.Zlimits, args.nTheta)

    if args.type == 'cone':
        if args.rlimits is None:
            print("ERROR: please specify minimum and maximum of cone radius",
                  file=sys.stderr)
            exit()

        flux = bin_cone(args.nZ, args.Zlimits, args.nTheta, args.rlimits)
        (x, y, z) = grid_cone(args.nZ, args.Zlimits, args.nTheta, args.rlimits)

    if args.type == 'sphere':
        (flux, R) = bin_sphere(args.nZ, args.Zlimits, args.nTheta)
        (x, y, z) = grid_sphere(R, args.nZ, args.Zlimits, args.nTheta)

    flux *= args.Pfactor
    plot_4D(flux, args.Fluxlimits, x, y, z)
