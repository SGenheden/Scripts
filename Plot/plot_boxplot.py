# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to plot boxplots
"""

import argparse

import numpy as np
import matplotlib.pylab as plt

from sgenlib import parsing

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description="Script to plot boxplots")
    argparser.add_argument('-f', '--files', nargs="+", help="the results")
    argparser.add_argument('-l', '--labels', help="the labels")
    argparser.add_argument('--ylabel',help="the y label")
    argparser.add_argument('--xlabels', nargs="+", help="the x labels", default=[])
    argparser.add_argument('-o', '--out', help="the output filename", default="box.png")
    argparser.add_argument('--rotate', action="store_true", help="rotate the boxplot", default=False)
    argparser.add_argument('--notches', action="store_true", help="add notches", default=False)
    argparser.add_argument('--squared', action="store_true",help="turn on squared plot",default=False)
    argparser.add_argument('--scaled', action="store_true",help="turn on scaled deviation",default=False)
    args = argparser.parse_args()

    error_list = []
    for filename in args.files :
        data = parsing.parse2ndarray(filename)
        error_list.append(np.abs(data[:,1] - data[:,0]))

    if args.labels is not None:
        print "\nOutliers (90 pecerntile):\n"
        with  open(args.labels, 'r') as f:
            labels = [line.strip() for line in f.readlines()]
        for i, errors in enumerate(error_list) :
            perc90= np.percentile(errors, 90)
            print ", ".join(["%s (%.1f)"%(label, error)
                for label,error in zip(labels, errors) if error >= perc90])
            print ""
            if args.scaled :
                error_list[i] = error_list[i]/error_list[i].max()

    if args.squared or args.correlation :
        f = plt.figure(figsize=(2.5,2.5))
        a = f.add_axes((0.2,0.2,0.75,0.75))
    else :
        f = plt.figure(figsize=(3.33,2.5))
        a = f.add_axes((0.15,0.2,0.75,0.75))
    a.boxplot(error_list,notch=args.notches, vert=(not args.rotate))
    if args.rotate :
        a.set_yticklabels(args.xlabels)
        if args.ylabel is not None : a.set_xlabel(args.ylabel,fontsize=8)
        if args.scaled : a.set_xlim([0,1.1])
    else :
        a.set_xticklabels(args.xlabels)
        if args.ylabel is not None : a.set_ylabel(args.ylabel,fontsize=8)
        if args.scaled : a.set_ylim([0,1.1])
    for tick in a.xaxis.get_major_ticks() :
        tick.label.set_fontsize(8)
    for tick in a.yaxis.get_major_ticks() :
        tick.label.set_fontsize(8)
    f.savefig(args.out,format="png",dpi=300)
