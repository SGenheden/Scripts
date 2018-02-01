# Author: Samuel Genheden, samuel.genheden@gmail.com

"""
Program to boostrap distribution
"""

import argparse
import math

import numpy as np
import matplotlib.pylab as plt

from sgenlib import parsing
from sgenlib import colors
from sgenlib import mol

if __name__ == '__main__' :

    parser = argparse.ArgumentParser(description="Boostrap a density file from g_density")
    parser.add_argument('-f', '--file', help="the density file from g_density")
    parser.add_argument('-d', '--dens', help="the density")
    parser.add_argument('--max',action="store_true", default=False, help="calculate the maxium position")
    parser.add_argument('--nboots', type=int, default=1000, help="the number of bootstraps")
    parser.add_argument('-o','--out',help="the output filename",default="boot.png")
    parser.add_argument('-a', '--area', type=float, help="membrane area")
    args = parser.parse_args()

    xvals, densities, ylabel = parsing.parse_densityfile(args.file)
    dx = (xvals[1] - xvals[0])
    lenz = xvals[-1]+dx-xvals[0]
    factor = len(xvals) / (args.area * lenz)
    midi = int(np.ceil(xvals.shape[0]/2.0))

    density = densities[args.dens] / factor
    if len(density) % 2 == 0 :
        dens = 0.5 * (density[:midi] + density[midi:][::-1])
    else :
        dens = 0.5 * (density[:midi] + density[midi-1:][::-1])
    xvals = xvals[:midi]
    if args.max :
        dens_std, max_std = mol.bootstrap_density(dens, xvals, nboots=args.nboots, calc_max=True)
        print "Maxium peak %.3f %.3f"%(xvals[np.argmax(dens)],max_std)
    else :
        dens_std = mol.bootstrap_density(dens, xvals, nboots=args.nboots)

    """f = plt.figure(figsize=(3.33,2.5))
    a = f.add_axes((0.15,0.2,0.75,0.75))
    a.plot(xvals, dens, "-",color=colors.color(0))
    plt.fill_between(xvals, dens-dens_std,dens+dens_std, facecolor='grey',linewidth=0,alpha=0.4,interpolate=True)
    f.savefig(args.out,format="png",dpi=300)"""
