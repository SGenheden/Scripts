# Author: Samuel Genheden, samuel.genheden@gmail.com

"""
Program to make 2D scatter plots of data
"""

import argparse

import numpy as np
import matplotlib.pylab as plt
import matplotlib.lines as mlines
import matplotlib.transforms as transforms

from sgenlib import colors
from sgenlib import parsing
from sgenlib import series

if __name__ == '__main__' :

    parser = argparse.ArgumentParser(description="Plot xy data")
    parser.add_argument('-f','--file',nargs="+",help="the data files")
    parser.add_argument('-l','--label',nargs="+",help="the labels")
    parser.add_argument('-o','--out',help="the output filename",default="xy.png")
    parser.add_argument('-s','--skip',type=float,help="how much to discard",default=0.0)
    parser.add_argument('--stride',type=int,help="how much to stride",)
    parser.add_argument('--xcol',type=int,help="the column of the x-data",default=0)
    parser.add_argument('--ycol',type=int,help="the column of the y-data",default=1)
    parser.add_argument('--multicol',nargs="+",type=int,help="the column of the y-data")
    parser.add_argument('--xfactor',type=float,default=1.0,help="multiplicator factor for x")
    parser.add_argument('--yfactor',type=float,default=1.0,help="multiplicator factor for y")
    parser.add_argument('--xlabel',help="the x label")
    parser.add_argument('--ylabel',help="the y label")
    parser.add_argument('--xmax',type=float,help="the maximum of the x-axis")
    parser.add_argument('--xmin',type=float,help="the minimum of the x-axis")
    parser.add_argument('--ymax',type=float,help="the maximum of the y-axis")
    parser.add_argument('--ymin',type=float,help="the minimum of the y-axis")
    parser.add_argument('--vertical',type=float,nargs="+",help="plot vertical lines",default=[])
    parser.add_argument('--vstyle', type=int, nargs="+",help="style of vertical lines")
    parser.add_argument('--equil',action="store_true",help="estimate equilibration",default=False)
    parser.add_argument('--nostyle',action="store_true",help="turn off different line styles",default=False)
    parser.add_argument('--scatter',action="store_true",help="turn on scatter line styles",default=False)
    parser.add_argument('--twin',action="store_true",help="turn on twin axis plotting",default=False)
    parser.add_argument('--errorbars',action="store_true",help="turn on error bar plotting",default=False)
    parser.add_argument('--squared', action="store_true",help="turn on squared plot",default=False)
    parser.add_argument('--correlation', action="store_true",help="turn on correlation line",default=False)
    parser.add_argument('--noyticks', action="store_true",help="turn on no y-ticks",default=False)
    parser.add_argument('--histogram', type=int, help="turn on plotting of histogram")
    parser.add_argument('--dotstyle', nargs="+",help="style for individual dots on the line")
    parser.add_argument('--xmarkers', help="put markers for the x-axis at the top of the graph")
    args = parser.parse_args()

    if len(args.file) > 2 :
        args.twin = False

    defaultstyle = "." if args.scatter else "-"

    if args.squared or args.correlation :
        f = plt.figure(figsize=(2.5,2.5))
        a = f.add_axes((0.19,0.19,0.74,0.74))
    else :
        f = plt.figure(figsize=(3.33,2.5))
        a = f.add_axes((0.19,0.19,0.74,0.74))
    plotax = [a]*len(args.file)
    if args.twin : plotax[1] = a.twinx()

    if args.ymax is not None and args.ymin is not None :
        for i, x in enumerate(args.vertical) :
            if args.vstyle is not None :
                stl = colors.style(args.vstyle[i])
            else :
                stl = colors.style(i+1)
            a.plot([x,x],[args.ymin,args.ymax], stl, color='k')

    ylim = [10E10, -10E10]
    xlim = [10E10, -10E10]

    lines = []
    for i, (filename,label, axis) in enumerate(zip(args.file,args.label, plotax)) :
        data = parsing.parse2ndarray(filename)
        if args.skip < 1 :
            n = int(np.floor(args.skip * data.shape[0]))
        else :
            n = args.skip
        data = data[n:,:]
        if args.multicol is not None :
            args.ycol = args.multicol[i]
        if data.shape[1] == 1 : args.ycol = 0
        y = data[n:,args.ycol]*args.yfactor
        if data.shape[1] > 1 :
            x = data[n:,args.xcol]*args.xfactor
        else:
            x = np.arange(1,len(y)+1)*args.xfactor
        if args.stride is not None :
            x = x[::args.stride]
            y = y[::args.stride]
        if args.histogram is not None:
            y, x = np.histogram(y, bins=args.histogram)
            x = 0.5 * (x[:-1] + x[1:])
        stl = defaultstyle if args.nostyle else colors.style(i)
        if args.nostyle and args.dotstyle :
            stl += args.dotstyle[min(i,len(args.dotstyle)-1)]
        ln = axis.plot(x,y, stl,color=colors.color(i),label=label)

        ylim[0] = min(ylim[0], y.min())
        ylim[1] = max(ylim[1], y.max())
        xlim[0] = min(xlim[0], y.min())
        xlim[1] = max(xlim[1], y.max())

        lines.append(ln[0])
        if args.errorbars:
            yerr = data[n:,args.ycol+1]*args.yfactor
            plt.fill_between(x,y-yerr,y+yerr, facecolor='grey',linewidth=0,alpha=0.4,interpolate=True)
        if args.equil :
            equili = series.find_equilibration(x,y, atleast=10)
            print "%s equilibrated at %.3f"%(label,x[equili])

    if len(args.file) == 2 and args.twin :
        a.legend(lines, args.label, loc='best', fancybox=True, framealpha=0.5,fontsize=8,labelspacing=0.20)
    elif len(args.file) > 1 :
        a.legend(loc='best', fancybox=True, framealpha=0.8,fontsize=8,labelspacing=0.20)
  #if len(args.file) > 1 : a.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0.0,fontsize=8)
    for ax in plotax :
        if args.noyticks :
            ax.set_yticks([])
        for tick in ax.xaxis.get_major_ticks() :
            tick.label.set_fontsize(8)
        for tick in ax.yaxis.get_major_ticks() :
            tick.label.set_fontsize(8)
            if args.twin :
                tick.label2.set_fontsize(8)
    if args.xlabel is not None : a.set_xlabel(args.xlabel,fontsize=8)
    if args.ylabel is not None : a.set_ylabel(args.ylabel,fontsize=8)

    if args.xmax is not None : a.set_xlim([0,args.xmax])
    if args.xmax is not None and args.xmin is not None : a.set_xlim([args.xmin,args.xmax])
    if args.ymax is not None and args.ymin is None : a.set_ylim([0,args.ymax])
    if args.ymax is not None and args.ymin is not None : a.set_ylim([args.ymin,args.ymax])
    if args.correlation or args.squared :
        if args.ymax is not None and args.ymin is not None :
            ylim[0] = args.ymin
            ylim[1] = args.ymax
        else:
            ylim[0] = min(ylim[0], xlim[0])
            ylim[1] = max(ylim[1], xlim[1])
        if args.correlation :
            plotax[0].plot([ylim[0],ylim[1]], [ylim[0],ylim[1]], '--k')
        a.set_ylim(ylim)
        a.set_xlim(ylim)
        a.set_aspect('equal')

    if args.xmarkers is not None :
        # Add a thin line to span the entire axis
        line = mlines.Line2D([0, 1], [1.02, 1.02], lw=1., color='k', clip_on=False, transform=a.transAxes)
        a.add_line(line)
        # Transformation to use x-coordinate in data coordinates, and
        # y-coordinate in axis coordinates
        trans = transforms.blended_transform_factory(a.transData, a.transAxes)
        # Add a bold line to marke a specific range of the x-axis
        with open(args.xmarkers, "r") as xfile :
            for line in xfile.readlines() :
                start, stop, lbl = line.strip().split()
                coord = map(float, [start, stop])
                # Bold line
                line = mlines.Line2D(coord, [1.02, 1.02], lw=3., color='k',
                                        clip_on=False, transform=trans)
                a.add_line(line)
                # Label
                a.text(np.mean(coord), 1.04, lbl, fontsize=6,
                        horizontalalignment='center', transform=trans)

    f.savefig(args.out,format="png",dpi=300)
