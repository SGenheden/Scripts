# Author: Samuel Genheden, samuel.genheden@gmail.com

"""
Program to calculate how many water molecules are leaking into the membrane
"""

import argparse
import math

import numpy as np

from sgenlib import parsing
from sgenlib import mol


def _count_water_inside(dens1, dens2, fi, li, fx, lx) :

    return sum(dens1[fi:li])

def _count_water_inside2(dens1, dens2, fi, li, fx, lx, **kwargs) :

    return sum(kwargs["watdens"][fi:li])

if __name__ == '__main__' :

    parser = argparse.ArgumentParser(description="Calculate water leakage")
    parser.add_argument('-f', '--file', help="the density file from g_density")
    parser.add_argument('-w', '--watdens', help="the water density")
    parser.add_argument('-l', '--lipdens', help="the lipid density")
    parser.add_argument('-i', '--lipdens2', help="the lipid density")
    parser.add_argument('-n', '--natom', type=int, help="the number of atoms in the lipid")
    parser.add_argument('-t', '--natom2', type=int, help="the number of atoms in the lipid")
    parser.add_argument('-a', '--area', type=float, help="membrane area")
    args = parser.parse_args()

    xvals, densities, ylabel = parsing.parse_densityfile(args.file)
    dx = (xvals[1] - xvals[0])
    lenz = xvals[-1]+dx-xvals[0]
    factor = len(xvals) / (args.area * lenz)

    densities[args.watdens] /= factor
    nwater = np.sum(densities[args.watdens])
    print "N water = %d "%(np.round(nwater))

    densities[args.lipdens] /= (factor * args.natom)
    nlipid = np.sum(densities[args.lipdens])
    print "N lipids = %d "%(np.round(nlipid))

    densities[args.lipdens2] /= (factor * args.natom2)
    nlipid = np.sum(densities[args.lipdens2])
    print "N lipids = %d "%(np.round(nlipid))

    midi = int(np.ceil(xvals.shape[0]/2.0))

    fi, li = mol.density_intercept(densities[args.watdens], densities[args.lipdens])
    (fx_std, lx_std), in_std = mol.bootstrap_density_intercept(densities[args.watdens],
            densities[args.lipdens], xvals, nboots=100, user_func=_count_water_inside)

    print "Water and lipid density crosses at: %.3f %.3f"%(xvals[fi], fx_std)
    print "Water and lipid density crosses at: %.3f %.3f"%(xvals[li], lx_std)
    print "\nNumber of leaked water: %d %d"%(_count_water_inside(
            densities[args.watdens], densities[args.lipdens],
            fi, li, xvals[fi], xvals[li]),
                in_std)

    fi, li = mol.density_intercept(densities[args.lipdens2], densities[args.lipdens])
    (fx_std, lx_std), in_std = mol.bootstrap_density_intercept(densities[args.lipdens2],
            densities[args.lipdens], xvals, nboots=100, user_func=_count_water_inside2,
            watdens=densities[args.watdens])

    print "\nLipid densities crosses at: %.3f %.3f"%(xvals[fi], fx_std)
    print "Lipid densities crosses at: %.3f %.3f"%(xvals[li], lx_std)
    print "\nNumber of leaked2 water: %d %d"%(_count_water_inside2(None, None,
            fi, li, xvals[fi], xvals[li], watdens=densities[args.watdens]),
                in_std)
