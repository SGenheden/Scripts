# Author: Samuel Genheden, samuel.genheden@gmail.com

"""
Program to make 2D scatter plots of densities computed with g_density
"""

import argparse
import re

import numpy as np
import matplotlib.pylab as plt

from sgenlib import colors
from sgenlib import parsing
from sgenlib import mol

if __name__ == '__main__' :

    parser = argparse.ArgumentParser(description="Plot density data")
    parser.add_argument('-f', '--file', nargs="+", help="the data files")
    parser.add_argument('-l', '--label', nargs="+", help="the labels")
    parser.add_argument('-d', '--densities', nargs="+", help="the densities to plot")
    parser.add_argument('-o', '--out', help="the output filename", default="densities.png")
    parser.add_argument('--scale', nargs="+",type=float, help="scale each density with this factor")
    parser.add_argument('--shiftdens', help="the density used to shift the x-axis")
    parser.add_argument('--half', action="store_true", help="only plot half the density", default=False)
    parser.add_argument('--xlabel',help="the x label")
    parser.add_argument('--ylabel',help="the y label")
    parser.add_argument('--noyticks', action="store_true",help="turn on no y-ticks",default=False)
    parser.add_argument('--nboots', type=int, default=0, help="the number of bootstraps")
    parser.add_argument('--max',action="store_true", default=False, help="calculate the maxium position")
    parser.add_argument('-a', '--area', type=float, nargs="+", help="membrane area")
    parser.add_argument('--dump', action="store_true",help="turn on dumping of data to standard out",default=False)
    args = parser.parse_args()

    f = plt.figure(figsize=(3.33,2.5))
    a = f.add_axes((0.15,0.2,0.75,0.75))

    plotthis = []
    for fi, filename in enumerate(args.file) :
        xvals, densities, ylabel = parsing.parse_densityfile(filename)
        midi = int(np.ceil(xvals.shape[0]/2.0))

        if args.area is not None :
            for d in densities :
                densities[d] /= mol.density_scaling(xvals, args.area[fi])

        if args.shiftdens is not None:
            maxi = np.argmax(densities[args.shiftdens][:midi])
            xvals = xvals - xvals[:midi][maxi]
        if args.half :
            xvals = xvals[:midi]

        for di, density in enumerate(args.densities):
            if len(args.file) == len(args.densities) :
                if di != fi : continue
                idx = di
            else :
                idx = di if len(args.file) == 1 else fi

            y = densities[density]
            if args.half :
                if len(y) % 2 == 0 :
                    y = 0.5 * (y[:midi] + y[midi:][::-1])
                else :
                    y = 0.5 * (y[:midi] + y[midi-1:][::-1])
#                y = y[:midi] # (y[:midi] + y[midi:][::-1]) * 0.5
            if args.scale is not None:
                y = y / args.scale[idx]

            if args.nboots > 0 :
                if args.max :
                    dens_std, max_std = mol.bootstrap_density(y, xvals, nboots=args.nboots, calc_max=True)
                else :
                    dens_std = mol.bootstrap_density(y, xvals, nboots=args.nboots)
                plt.fill_between(xvals, y-dens_std,y+dens_std, facecolor='grey',linewidth=0,alpha=0.4,interpolate=True)
            a.plot(xvals , y, "-", color=colors.color(idx),label=args.label[idx])

            if args.dump :
                for xi, yi in zip(xvals, y):
                    print "%.3E %.3E"%(xi, yi)
                print "---"

            if args.max :
                dmax = xvals[np.argmax(y)]
                if args.nboots > 0 :
                    print "Maxium peak at\t%.3f\t%.3f\twhich is %.3f\t%.3f"%(dmax, max_std, y.max(), dens_std[np.argmax(y)])
                else :
                    print "Maxium peak at\t%.3f\twhich is %.3f"%(dmax, y.max())


    if len(args.file) > 1 or len(args.densities) > 1 :
        a.legend(loc='best', fancybox=True, framealpha=0.5,fontsize=8,labelspacing=0.20)

    for tick in a.xaxis.get_major_ticks() :
        tick.label.set_fontsize(8)
    for tick in a.yaxis.get_major_ticks() :
        tick.label.set_fontsize(8)
    if args.xlabel is None :
        if args.shiftdens is not None:
            a.set_xlabel("Distance from %s [nm]"%args.shiftdens,fontsize=8)
        else:
            a.set_xlabel("z [nm]",fontsize=8)
    else :
        a.set_xlabel(args.xlabel,fontsize=8)
    if args.ylabel is None :
        a.set_ylabel(ylabel,fontsize=8)
    else:
        a.set_ylabel(args.ylabel,fontsize=8)
    if args.noyticks :
        a.set_yticks([])
    f.savefig(args.out,format="png",dpi=300)
