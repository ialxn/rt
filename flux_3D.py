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
import pandas as pd
import sys
import matplotlib.pyplot as plt


def bin(P_per_area, nZ, nTheta, dZ, dTheta, H):
    """Bin x,y,z data in data frame df

    Using numpy.histogram2d() does not provide any speed benefit but we
    loose the possibility to track the number of missed tuples.

    Args:
        df: DataFrame holding carthesian coordinates
        P_per_area: Flux density corresponding to one absorbed ray
        nZ (int): Number of bins in z-direction
        nTheta: Number of bins in Theta-direction
        dZ: bin width in z-direction
        dTheta: bin width in Theta direction
        H: z_min, z_max

    Returns:
        binned (nZ * nTheta+1): binned flux data
    """
    df = pd.read_csv(sys.stdin, delimiter='\s+',
                     header=None, names=['x', 'y', 'z'], usecols=(0, 1, 2))

    binned = np.zeros((nZ, nTheta + 1))
    s = 0
    for (X, Y, Z) in zip(df.x, df.y, df.z):
        theta = np.pi + np.arctan2(Y, X)   # theta 0..2pi to avoid negative j
        i = int(np.trunc((Z - H[0]) / dZ))
        j = int(np.trunc(theta / dTheta))
        # precautions for round off errors
        if i >= nZ:
            s += 1
            continue
        if j > nTheta:
            s += 1
            continue
        binned[i, j] += P_per_area[i]
    # DataFrame size correct for header line
    print("{0} of {1} tuples skipped.".format(s, df.x.size - 1))

    return binned


def bin_cone(n_z, H, n_theta, r):
    """Bin flux data on cone

    Args:
        n_z: Number of bins in z-direction
        H: Minimum and Maximum z values
        n_theta: number of bins in theta direction
        r: Radius at base and at top of cone

    Returns:
        binned: Binned flux data on a nZ * n_theta+1 grid
        R_xy: Radius of cone as function of z-coordinate
    """

    L = H[1] - H[0]
    Dz = L / n_z    # intervall in z direction
    Dtheta = 2 * np.pi / (n_theta+1)
    R = r[0] - r[1]
    dl = np.sqrt(L**2 + R**2)   # length of line on surface of cone
    R_xy = np.linspace(r[0] + Dz * R / (2 * L),
                       r[0] + Dz * R / (2 * L) - n_z * Dz * R / L,
                       n_z)

    area = np.abs(R_xy * Dtheta * dl)
    binned = bin(1/area, n_z, n_theta, Dz, Dtheta, H)

    return binned, R_xy


def bin_cylinder(n_z, H, n_theta):
    """Bin flux data on cylinder

    Args:
        n_z: Number of bins in z-direction
        H: Minimum and Maximum z values
        n_theta: number of bins in theta direction

    Returns:
        binned: Binned flux data on a nZ * n_theta+1 grid
        R: Radius of cylinder as function of z-coordinate
    """
    l = sys.stdin.readline().strip().split()
    x = float(l[0])
    y = float(l[1])
    r = np.sqrt(x**2 + y**2)

    R = np.asarray([r for i in range(n_z)])
    Dz = (H[1] - H[0]) / n_z
    Dtheta = 2 * np.pi / (n_theta + 1)
    area = np.abs(R * Dtheta * Dz)

    binned = bin(1/area, n_z, n_theta, Dz, Dtheta, H)

    return binned, R


def bin_paraboloid(n_z, H, n_theta, foc):
    """Bin flux data on paraboloid

    Args:
        n_z: Number of bins in z-direction
        H: Minimum and Maximum z values
        n_theta: number of bins in theta direction
        foc: focal length of paraboloid

    Returns:
        binned: Binned flux data on a nZ * n_theta+1 grid
        R_xy: Radius of paraboloid as function of z-coordinate
    """
    Dz = (H[1] - H[0]) / n_z
    Dtheta = 2 * np.pi / (n_theta + 1)

    h = np.linspace(H[0], H[1], n_z + 1)
    R_xy = np.sqrt(4 * foc * h[:n_z])

    ap = 1.0 / (4 * foc)
    A0 = np.sqrt(4 * h[:n_z] / ap + 1 / ap**2) * (1 + 4 * ap * h[:n_z])
    A1 = np.sqrt(4 * h[1:] / ap + 1 / ap**2) * (1 + 4 * ap * h[1:])
    area = np.pi / (6 * ap) * (A1 - A0) / n_theta

    binned = bin(1/area, n_z, n_theta, Dz, Dtheta, H)

    return binned, R_xy


def bin_sphere(n_z, H, n_theta):
    """
    Read output of 'get_flux' from stdin.
    Perform binning of the data on the grid 'n_z' by 'n_theta+1'.

        H: (min_height, max_height)

    return the binned data (flux density)
    """
    l = sys.stdin.readline().strip().split()
    x = float(l[0])
    y = float(l[1])
    z = float(l[2])
    R = np.sqrt(x**2 + y**2 + z**2)

    Dz = (H[1] - H[0]) / n_z
    Dtheta = 2 * np.pi / (n_theta + 1)
    h = np.linspace(H[0], H[1], n_z + 1) / R
    phi1 = np.arccos(h[:n_z])
    Dphi = np.arccos(h[1:]) - phi1
    R_xy = R * np.sin(phi1)

    area = np.abs(R**2 * Dphi * Dtheta)
    binned = bin(1/area, n_z, n_theta, Dz, Dtheta, H)

    return binned, R_xy


def R2xy(R, theta):
    """Convert polar to carthesian coordinates

    Args:
        R: Radius
        theta: Angle

        Arguments can be scalars or vectors

    Returns:
        x, y: carthesian coordinates
    """
    x = R * np.sin(theta)
    y = R * np.cos(theta)

    return x, y


def grid(R, nZ, H, nTheta):
    """Calculate carthesian coordinates on grid (nZ * nTheta+1)

    Calculates carthesian coordinates of rotationally symmetrical body
    (rotated around z-axis) defined by R(z) between H[0] and H[1].

    Args:
        R (float[nZ]): Radius as function of z
        nZ (int): Number of grid points in z direction
        H (float,float): minimum and maximum z values
        nTheta (int): Number of grid points in direction theta

    Returns:
        x, y, z: Carthesian coordinates on grid nz x ntheta+1)
    """
    z = np.linspace(H[0], H[1], nZ)
    theta = np.linspace(-np.pi, np.pi, nTheta + 1)
    Ri = np.transpose(np.asarray([R for i in range(nTheta + 1)]))

    thetai, zi = np.meshgrid(theta, z)
    xi, yi = R2xy(Ri, thetai)

    return xi, yi, zi


def axisEqual3D(ax):
    """Sets aspect ratio between x,y,z to 1.0

    Args:
        ax: Axes instance

    Returns:
        None
    """
    exts = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    sz = exts[:, 1] - exts[:, 0]
    centers = np.mean(exts, axis=1)
    maxsize = max(abs(sz))
    r = maxsize/2
    for ctr, dim in zip(centers, 'xyz'):
        getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)


def plot_4D(flux, Limits, x, y, z):
    """Plots flux distribution on 3D body

    Plots (colored patches using Limits as minimum and maximum values) flux
    distribution on body defined by the carthesian coordinates (x,y,z).
    flux, x, y, z have to be on identical grid (nZ * nTheta+1). If Limits is
    None autoscale is applied.

    Args:
        flux: flux values (on grid nZ * nTheta+1)
        Limits (float,float): Minimum and maximum flux values to be used to
                              set colorscale. If None autoscale is performed
        x, y, z,: Carthesian coordinates (on grid nZ * nTheta+1)

    Returns:
        plt: Plot
    """
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm, colors

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

    supported = ','.join(plt.figure().canvas.get_supported_filetypes())

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
                        choices=['cone', 'cylinder', 'paraboloid', 'sphere'],
                        required=True,
                        help='Type of body')
    parser.add_argument('-r', '--rlimits', metavar='R', nargs=2, type=float,
                        help='Minimum and maximum value of cone radius')
    parser.add_argument('-f', '--foc', metavar='f', type=float,
                        help='Focal length of paraboloid')
    parser.add_argument('-P', '--Pfactor', metavar='P', type=float,
                        default=1.0,
                        help='Scaling factor for flux "P_factor"')
    parser.add_argument('-o', '--output', metavar='FILE', type=str,
                        help='Output to file. Supported formats: ' + supported)

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

    elif args.type == 'paraboloid': ## PARABOLOID ####
        if args.foc is None:
            print("ERROR: please specify focal length of paraboloid",
                  file=sys.stderr)
            exit()
        (flux, R) = bin_paraboloid(args.nZ, args.Zlimits, args.nTheta, args.foc)
        (x, y, z) = grid(R, args.nZ, args.Zlimits, args.nTheta)

    elif args.type == 'sphere': ###### SPHERE ########
        (flux, R) = bin_sphere(args.nZ, args.Zlimits, args.nTheta)
        (x, y, z) = grid(R, args.nZ, args.Zlimits, args.nTheta)

    flux *= args.Pfactor
    plt = plot_4D(flux, args.Fluxlimits, x, y, z)

    if args.output is None:
        plt.show()
    else:
        try:
            plt.savefig(args.output)
        except ValueError:
            print('Cannot save figure ({})'.format(args.output))
            print('Supported formats: {}'.format(supported))
            plt.show()
