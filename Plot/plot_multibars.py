# Author: Samuel Genheden, samuel.genheden@gmail.com

"""
Program to make bar plots of data

DOES NOT work great

"""

import argparse

import numpy as np
import matplotlib.pylab as plt

from sgenlib import colors
from sgenlib import parsing

def _plot_single(data, errors, headers, barlabels, ylabel, fig):

    a = fig.add_axes((0.2,0.15,0.75,0.75))

    ngroups = data.shape[0]
    nitems = len(headers)
    print ngroups, nitems

    gleft = 0.8*(np.arange(nitems)+1.0)
    width = gleft[-1] - gleft[0]
    left =  np.asarray([gleft+(width + 1.6)*i for i in range(ngroups)]).reshape(ngroups*nitems)
    gcolor = [colors.color(i) for i in range(nitems)]
    color = []
    labels = []
    for i in range(ngroups):
        color.extend(gcolor)
        labels.extend("")
    a.bar(left, data.reshape(ngroups*nitems), yerr=errors.reshape(ngroups*nitems),
            width=0.8, color=color, error_kw=dict(ecolor='k'), label=headers)
    a.set_xticks(left[::2]+width)
    a.set_xticklabels(barlabels)
    #a.legend()
    h, l = a.get_legend_handles_labels()
    a.legend(h[0], headers, loc='best', fancybox=True, framealpha=0.5,fontsize=8,labelspacing=0.20)

    for tick in a.xaxis.get_major_ticks() :
        tick.label.set_fontsize(8)
    for tick in a.yaxis.get_major_ticks() :
        tick.label.set_fontsize(8)
    a.set_ylim([data.min()-3.0,data.max()+3.0])
    #for i,h  in enumerate(barlabels):
    #    a.text(left[nitems*i]+0.4,data.min()-6.0,h, size=8)
    a.set_ylabel(ylabel, fontsize=8)
    a.set_xlabel("Docking pose", fontsize=8)

def _plot_multi(data, errors, headers, barlabels, nrow, ncol, fig):

    for i, (idata, ierrors, header) in enumerate(zip(data.T, errors.T, headers), 1):
        a = fig.add_subplot(nrow, ncol, i)
        left = (np.arange(data.shape[0])+1)*2-0.4
        barcolors = [colors.color(i) for i in range(data.shape[0])]
        a.bar(left, idata, yerr=ierrors, color=barcolors, error_kw=dict(ecolor='k'))
        a.set_xticks(left+0.4)
        a.set_xticklabels(barlabels, rotation=45)
        a.set_ylabel(header, fontsize=8)
        for tick in a.xaxis.get_major_ticks() :
            tick.label.set_fontsize(8)
        for tick in a.yaxis.get_major_ticks() :
            tick.label.set_fontsize(8)

if __name__ == '__main__' :

    parser = argparse.ArgumentParser(description="Plot bar plots")
    parser.add_argument('-f','--file',help="the data for the bar plots")
    parser.add_argument('-y','--ylabel',help="the ylabel", default="")
    parser.add_argument('--single', action="store_true", help="plot all in one sub-plot", default=False)
    parser.add_argument('-o','--out',help="the output filename",default="bar.png")
    args = parser.parse_args()

    with open(args.file,"r") as f :
        lines = f.readlines()
    headers = lines[0].strip().split("\t")[::2]
    barlabels = []
    data = []
    for line in lines[1:]:
        cols = line.strip().split("\t")
        barlabels.append(cols[0])
        data.append(cols[1:])
    data = np.array(data, dtype=float)
    errors = data[:,1::2]
    data = data[:,::2]
    print data, errors, headers, barlabels

    if args.single :
        ncol = 1
        nrow = 1
    else :
        if len(headers) == 2 :
            ncol = 1
            nrow = 2
        else:
            ncol = 2
            nrow = int(np.ceil(len(headers)/2.0))
    widths = [0,3.25,7]
    f = plt.figure(0, figsize=(widths[ncol], 2.5*nrow))

    if args.single :
        _plot_single(data, errors, headers, barlabels, args.ylabel, f)
    else :
        _plot_multi(data, errors, headers, barlabels, nrow, ncol, f)
    #f.tight_layout()
    f.savefig(args.out, format="png", dpi=300)
