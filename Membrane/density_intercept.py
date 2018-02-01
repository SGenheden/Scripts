# Author: Samuel Genheden, samuel.genheden@gmail.com

"""
Program to compute where two densities cross
"""

import argparse
import math

import numpy as np

from sgenlib import parsing

if __name__ == '__main__' :

    parser = argparse.ArgumentParser(description="Calculate density intercept")
    parser.add_argument('-f', '--file', help="the density file from g_density")
    parser.add_argument('-s', '--scaling', type=float, nargs=2, help="the scaling factor for the densities", default=[1.0,1.0])
    parser.add_argument('-d', '--dens', nargs=2, help="the densities")
    args = parser.parse_args()

    xvals, densities, ylabel = parsing.parse_densityfile(args.file)
    dens1 = densities[args.dens[0]] / args.scaling[0]
    dens2 = densities[args.dens[1]] / args.scaling[1]

    midi = int(np.ceil(xvals.shape[0]/2.0))

    i = 0
    while dens1[i] == 0.0 or dens2[i] < dens1[i] :
        i += 1
    print "%s and %s crosses at: %.3f"%(args.dens[0], args.dens[1], xvals[i])

    i = len(xvals) - 1
    while  dens1[i] == 0.0 or dens2[i] < dens1[i] :
        i -= 1
    print "%s and %s crosses at: %.3f"%(args.dens[0], args.dens[1], xvals[i])
