# Author: Samuel Genheden, samuel.genheden@gmail.com

"""
Program to make 2D scatter plots of densities computed with g_density
"""

import argparse
import re

import numpy as np
import matplotlib.pylab as plt

from sgenlib import colors

def _parse_xvgfile(filename) :

    lines = []
    with open(filename, 'r') as f:
        lines = f.readlines()

    i = 0
    while lines[i].startswith("#") : i = i + 1

    legends = []
    while lines[i].startswith('@'):
        if lines[i].find("yaxis  label") > 0 :
            cols = lines[i].split('"')
            ylabel = cols[-2]
        elif re.search(" s\d+ legend",lines[i]):
            cols = lines[i].split('"')
            legends.append(cols[-2])
        i = i + 1

    data = []
    while i < len(lines):
        data.append(lines[i].strip().split())
        i = i + 1
    data = np.array(data, dtype=float)

    densities = {l : col for l, col in zip(legends, data[:,1:].T)}
    return data[:,0], densities, ylabel

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
    args = parser.parse_args()

    f = plt.figure(figsize=(3.33,2.5))
    a = f.add_axes((0.15,0.2,0.75,0.75))

    plotthis = []
    for fi, filename in enumerate(args.file) :
        xvals, densities, ylabel = _parse_xvgfile(filename)
        midi = int(np.ceil(xvals.shape[0]/2.0))
        if args.shiftdens is not None:
            maxi = np.argmax(densities[args.shiftdens][:midi])
            xvals = xvals - xvals[:midi][maxi]
        if args.half :
            xvals = xvals[:midi]
        for di, density in enumerate(args.densities):
            idx = di if len(args.file) == 1 else fi
            y = densities[density]
            if args.half :
                y = y[:midi]
            if args.scale is not None:
                y = y / args.scale[idx]
            a.plot(xvals , y, "-", color=colors.color(idx),label=args.label[idx])

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
