# Author: Samuel Genheden, samuel.genheden@gmail.com

"""
Program to make bar plots of data
"""

import argparse

import numpy as np
import matplotlib.pylab as plt

from sgenlib import colors
from sgenlib import parsing

if __name__ == '__main__' :

    parser = argparse.ArgumentParser(description="Plot bar plots")
    parser.add_argument('-f','--file',help="the data for the bar plots")
    parser.add_argument('-o','--out',help="the output filename",default="bar.png")
    args = parser.parse_args()

    with open(args.file,"r") as f :
        lines = f.readlines()
    headers = lines[0].strip().split("\t")[::2]
    barlabels = []
    data = []
    for line in lines[1:]:
        cols = line.strip().split()
        barlabels.append(cols[0])
        data.append(cols[1:])
    data = np.array(data, dtype=float)
    errors = data[:,1::2]
    data = data[:,::2]

    if len(headers) == 2 :
        ncol = 1
        nrow = 2
    else:
        ncol = 2
        nrow = int(np.ceil(len(headers)/2.0))
    widths = [0,3.25,7]

    f = plt.figure(0, figsize=(widths[ncol], 2.5*nrow))

    for i, (idata, ierrors, header) in enumerate(zip(data.T, errors.T, headers), 1):
        a = f.add_subplot(nrow, ncol, i)
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
    f.tight_layout()
    f.savefig(args.out, format="png", dpi=300)
