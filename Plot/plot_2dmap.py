# Author: Samuel Genheden, samuel.genheden@gmail.com

"""
Program to plot 2D maps. The maps are named after the input files but
with a "png"-extension.

Example:
    plot_2dmap.py -f grid_dhh_low.dat grid_dhh_upp.dat
"""

import argparse
import os

import numpy as np
import matplotlib.pylab as plt

from sgenlib import parsing

def _draw_2dmap(mat, cutoff, filename) :

    f = plt.figure()
    sel = mat > cutoff
    minval = np.min(mat[sel])
    im = f.gca().imshow(mat, vmin=minval, cmap=plt.cm.RdYlBu,origin="lower")
    excluded_im = np.ones([mat.shape[0],mat.shape[1],4])
    excluded_im[sel,3] = 0.0
    f.gca().imshow(excluded_im, origin="lower")
    f.colorbar(im)
    f.savefig(filename,format="png",dpi=300)

if __name__ == '__main__' :

    parser = argparse.ArgumentParser(description="Plotting 2D maps")
    parser.add_argument('-f', '--file',nargs="+",help="the input file")
    parser.add_argument('--cutoff', type=float, default=0.0, help="the cutoff for filtering")
    args = parser.parse_args()

    for filename in args.file :
        mat = parsing.parse2ndarray(filename)
        _draw_2dmap(mat, args.cutoff, os.path.splitext(filename)[0]+".png")
