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


def bin(x, y, z, P_per_area, nZ, nTheta, dZ, dThetai, H):
    binned = np.zeros((nZ, nTheta + 1))
    for (X, Y, Z) in zip(x, y, z):
        theta = np.pi + np.arctan2(Y, X)   # theta 0..2pi to avoid negative j
        i = int(np.trunc((Z - H[0]) / dZ))
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

    L = H[1] - H[0]
    Dz = L / n_z    # intervall in z direction
    Dtheta = 2 * np.pi / (n_theta+1)
    R = r[0] - r[1]
    dl = np.sqrt(L**2 + R**2)   # length of line on surface of cone
    R_xy = np.linspace(r[0] + Dz * R / (2 * L),
                       r[0] + Dz * R / (2 * L) - n_z * Dz * R / L,
                       n_z)

    area = np.abs(R_xy * Dtheta * dl)
    binned = bin(x, y, z, 1/area, n_z, n_theta, Dz, Dtheta, H)

    return binned, R_xy


def bin_cylinder(n_z, H, n_theta):
    """
    Read output of 'get_flux' from stdin.
    Perform binning of the data on the grid 'n_z' by 'n_theta+1'.

        H: (0, height_of_cylinder)

    return the binned data (flux density) and
    the scalar 'R' (radius of cylinder)
    """
    (x, y, z) = np.loadtxt(sys.stdin, usecols=(0,1,2), unpack=True)

    Dz = (H[1] - H[0]) / n_z
    Dtheta = 2 * np.pi / (n_theta + 1)
    r = np.sqrt(x[0]**2 + y[0]**2)
    R = np.asarray([r for i in range(n_z)])

    area = np.abs(R * Dtheta * Dz)
    binned = bin(x, y, z, 1/area, n_z, n_theta, Dz, Dtheta, H)

    return binned, R


def bin_sphere(n_z, H, n_theta):
    """
    Read output of 'get_flux' from stdin.
    Perform binning of the data on the grid 'n_z' by 'n_theta+1'.

        H: (min_height, max_height)

    return the binned data (flux density)
    """
    (x, y, z) = np.loadtxt(sys.stdin, usecols=(0,1,2), unpack=True)

    Dz = (H[1] - H[0]) / n_z
    Dtheta = 2 * np.pi / (n_theta + 1)
    R = np.sqrt(x[0]**2 + y[0]**2 + z[0]**2)
    h = np.linspace(H[0], H[1], n_z + 1) / R
    phi1 = np.arccos(h[:n_z])
    Dphi = np.arccos(h[1:]) - phi1
    R_xy = R * np.sin(phi1)

    area = np.abs(R**2 * Dphi * Dtheta)
    binned = bin(x, y, z, 1/area, n_z, n_theta, Dz, Dtheta, H)

    return binned, R_xy


def R2xy(R, theta):
    """
    converts 'R', 'theta' to 'x','y'
    'R' can be scalar or vector of the same size as axis0 of 'theta'
    """
    x = R * np.sin(theta)
    y = R * np.cos(theta)

    return x, y


def grid(R, nZ, H, nTheta):
    z = np.linspace(H[0], H[1], nZ)
    theta = np.linspace(0, 2 * np.pi, nTheta + 1)
    Ri = np.transpose(np.asarray([R for i in range(nTheta + 1)]))

    thetai, zi = np.meshgrid(theta, z)
    xi, yi = R2xy(Ri, thetai)

    return xi, yi, zi


def axisEqual3D(ax):
    exts = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    sz = exts[:, 1] - exts[:, 0]
    centers = np.mean(exts, axis=1)
    maxsize = max(abs(sz))
    r = maxsize/2
    for ctr, dim in zip(centers, 'xyz'):
        getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)


def plot_4D(flux, Limits, x, y, z):
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm, colors
    import matplotlib.pyplot as plt

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
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

    axisEqual3D(ax)

    cbar = fig.colorbar(s_m)
    cbar.set_label('Flux')

    return plt

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

    if args.type == 'cylinder': ###### CYLINDER ######
        (flux, R) = bin_cylinder(args.nZ, args.Zlimits, args.nTheta)
        (x, y, z) = grid(R, args.nZ, args.Zlimits, args.nTheta)

    elif args.type == 'cone': ######## CONE ##########
        if args.rlimits is None:
            print("ERROR: please specify minimum and maximum of cone radius",
                  file=sys.stderr)
            exit()

        (flux, R) = bin_cone(args.nZ, args.Zlimits, args.nTheta, args.rlimits)
        (x, y, z) = grid(R, args.nZ, args.Zlimits, args.nTheta)

    elif args.type == 'sphere': ###### SPHERE ########
        (flux, R) = bin_sphere(args.nZ, args.Zlimits, args.nTheta)
        (x, y, z) = grid(R, args.nZ, args.Zlimits, args.nTheta)

    flux *= args.Pfactor
    plt = plot_4D(flux, args.Fluxlimits, x, y, z)

    plt.show()
