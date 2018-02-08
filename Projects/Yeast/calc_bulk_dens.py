# Author: Samuel Genheden, samuel.genheden@gmail.com

"""
Program to compute bulk density of solutes from densities
"""

import argparse
import math
import re

import numpy as np

from sgenlib import parsing
from sgenlib import mol

def _count_solute_inside(dens1, dens2, fi, li, fx, lx, **kwargs) :

    return sum(kwargs["soldens"][fi:li])

def _membrane_vol(dens1, dens2, fi, li, fx, lx, **kwargs) :

    return (lx - fx)*kwargs["area"]

if __name__ == '__main__' :

    parser = argparse.ArgumentParser(description="Calculate bulk density")
    parser.add_argument('-f', '--file', help="the density file from g_density")
    parser.add_argument('-s', '--soldens', help="the solute density")
    parser.add_argument('-w', '--watdens', help="the water density")
    parser.add_argument('-l', '--lipdens', help="the lipid density")
    parser.add_argument('-n', '--natom', type=int, help="the number of atoms in the lipid")
    parser.add_argument('-a', '--area', type=float, help="membrane area")
    args = parser.parse_args()

    xvals, densities, ylabel = parsing.parse_densityfile(args.file)
    dx = (xvals[1] - xvals[0])
    lenz = xvals[-1]+dx-xvals[0]
    factor = len(xvals) / (args.area * lenz)

    densities[args.soldens] /= factor
    nsolute = np.sum(densities[args.soldens])
    print "N solutes = %d "%(np.round(nsolute))

    densities[args.watdens] /= factor
    nwater = np.sum(densities[args.watdens])
    print "N water = %d "%(np.round(nwater))

    densities[args.lipdens] /= (factor * args.natom)
    nlipid = np.sum(densities[args.lipdens])
    print "N lipids = %d "%(np.round(nlipid))

    midi = int(np.ceil(xvals.shape[0]/2.0))

    """fi = 0
    while  densities[args.lipdens][fi] / nlipid < densities[args.watdens][fi] / nwater :
        fi += 1
    print "Water and lipid density crosses at: %.3f"%(xvals[fi])
    nfirst = np.sum(densities[args.soldens][:fi])
    nfirst_wat = np.sum(densities[args.watdens][:fi])
    print "N solutes = %d in %d waters"%(np.round(nfirst), np.round(nfirst_wat))

    li = len(xvals) - 1
    while  densities[args.lipdens][li] / nlipid < densities[args.watdens][li] / nwater :
        li -= 1
    print "Water and lipid density crosses at: %.3f"%(xvals[li])
    nsecond = np.sum(densities[args.soldens][li:])
    nsecond_wat = np.sum(densities[args.watdens][li:])
    print "N solutes = %d in %d waters"%(np.round(nsecond), np.round(nsecond_wat))"""

    fi, li = mol.density_intercept(densities[args.watdens], densities[args.lipdens])
    (fx_std, lx_std), nsol_std = mol.bootstrap_density_intercept(densities[args.watdens],
            densities[args.lipdens], xvals, nboots=100, user_func=_count_solute_inside,
                soldens=densities[args.soldens])

    nfirst = np.sum(densities[args.soldens][:fi])
    nfirst_wat = np.sum(densities[args.watdens][:fi])
    nsecond = np.sum(densities[args.soldens][li:])
    nsecond_wat = np.sum(densities[args.watdens][li:])
    print "Water and lipid density crosses at: %.3f %.3f"%(xvals[fi], fx_std)
    print "Water and lipid density crosses at: %.3f %.3f"%(xvals[li], lx_std)
    print "N solutes = %d in %d waters"%(np.round(nfirst), np.round(nfirst_wat))
    print "N solutes = %d in %d waters"%(np.round(nsecond), np.round(nsecond_wat))

    # Estimate how much of the simulation volume is water
    (fx_std, lx_std), vol_std = mol.bootstrap_density_intercept(densities[args.watdens],
            densities[args.lipdens], xvals, nboots=100, user_func=_membrane_vol,
                area=args.area)
    mem_vol = _membrane_vol(None, None, fi, li, xvals[fi], xvals[li], area=args.area)
    all_vol = lenz * args.area
    wat_vol = all_vol - mem_vol
    first_vol = xvals[fi] * args.area
    last_vol  = (lenz - xvals[li]) * args.area
    all_vol = (xvals[-1] + xvals[1] - xvals[0]) * args.area
    print "Water volume %.3f"%wat_vol
    print "Total volume %.3f"%all_vol
    print "Water content %.3f"%(wat_vol/all_vol)
    print "Membrane volume %.3f"%mem_vol

    # Now for bulk concentration calculation
    avogrado = 6.022 * math.pow(10.0,23.0)
    if args.soldens[:3] == "eth" :
        molmass = 46.07 # g / mol
        dens = 789.0 # g / L
    elif args.soldens[:3] == "but" :
        molmass = 74.12 # g / mol
        dens = 810.0 # g / L
    conc_factor = molmass / (avogrado * math.pow(10.0,-24.0))

    nsolute_bulk = float(nfirst + nsecond)
    konc = nsolute_bulk / wat_vol
    konc_std = np.abs(konc)*np.sqrt(math.pow(nsol_std/nsolute_bulk,2)+math.pow(vol_std/wat_vol,2))
    print "Bulk concentration = %.3f %.3f"%(conc_factor*konc, conc_factor*konc_std)

    # Lipid concentration
    nsolute_lip = float(nsolute - nfirst - nsecond)
    konc = nsolute_lip / mem_vol
    konc_std = np.abs(konc)*np.sqrt(math.pow(nsol_std/nsolute_lip,2)+math.pow(vol_std/mem_vol,2))
    print "Lipid concentration = %.3f %.3f"%(conc_factor*konc, conc_factor*konc_std)

    print "Fraction in membrane = %.3f %.3f"%(nsolute_lip/nsolute,  nsol_std/nsolute)
