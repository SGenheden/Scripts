# Author: Samuel Genheden, samuel.genheden@gmail.com

"""
Program to plot 2D maps. The maps are named after the input files but
with a "png"-extension.

Example:
    plot_2dmap.py -f grid_dhh_low.dat grid_dhh_upp.dat
    plot_2dmap.py -f grid_dhh_low.dat grid_dhh_upp.dat --ylabel "Y [nm]" --xlabel "X [nm]" --extent 0 10 0 10 --cblabel "Dhh [nm]" --max 5.0
    plot_2dmap.py  -f grid_dhh_low.dat grid_dhh_upp.dat --ylabel "Y [nm]" --xlabel "X [nm]" --extent 0 10 0 10 --cblabel "Dhh [nm]" --multiplot -o grid_dhh.png
"""

import argparse
import string
import os

import numpy as np
import matplotlib.pylab as plt
import matplotlib.image as mimage

from sgenlib import parsing

def _draw_mat(mat, cutoff, maxval, xlabel, ylabel, extent, axis) :

    sel = mat > cutoff
    minval = np.min(mat[sel])
    im = axis.imshow(mat, vmin=minval, vmax=maxval, cmap=plt.cm.RdYlBu,origin="lower",extent=extent)
    excluded_im = np.ones([mat.shape[0],mat.shape[1],4])
    excluded_im[sel,3] = 0.0
    axis.imshow(excluded_im, origin="lower",extent=extent)
    for tick in axis.xaxis.get_major_ticks() :
        tick.label.set_fontsize(8)
    for tick in axis.yaxis.get_major_ticks() :
        tick.label.set_fontsize(8)

    if xlabel is not None : axis.set_xlabel(xlabel,fontsize=8)
    if ylabel is not None : axis.set_ylabel(ylabel,fontsize=8)

    return im

def _draw_2dmap(mat, cutoff, maxval, xlabel, ylabel, extent, cblabel, nobar, filename) :

    if nobar :
        f = plt.figure(figsize=(2.5,2.5))
        a = f.add_axes((0.2,0.15,0.75,0.75))
    else :
        f = plt.figure(figsize=(3.33,2.5))
        a = f.add_axes((0.15,0.15,0.75,0.75))

    im = _draw_mat(mat, cutoff, maxval, xlabel, ylabel, extent, a)

    if not nobar :
        cb = f.colorbar(im, pad=0.1)
        cb.ax.tick_params(labelsize=8)
        if cblabel is not None : cb.ax.set_title(cblabel, size=8)

    f.savefig(filename,format="png",dpi=300)

    return im

def _draw_multi_2dmap(mats, cutoff, maxval, ylabel, xlabel, extent, cblabel, filename) :

    if len(mats) > 4 or len(mats) == 1 :
        raise Exception("Multiplot is only possible with 2 - 4 matrices")

    cols = 2
    rows = 2 if len(mats) > 2 else 1
    f = plt.figure(figsize=(3.33*cols,2.5*rows))
    f.subplots_adjust(wspace=0.0,hspace=0,left=0.00,right=1.0,top=1.0,bottom=0.00)

    labels = [c+")" for c in string.ascii_uppercase[:len(args.file)]]
    props = dict(facecolor='white',edgecolor='white', alpha=1.0)
    if maxval is None :
        maxval = np.max(np.asarray(mats))
    axes = []
    for i, (mat, label) in enumerate(zip(mats, labels), 1) :
        a = f.add_subplot(rows, cols,i)
        axes.append(a)
        raw_im = _draw_2dmap(mat, cutoff, maxval, xlabel, ylabel, extent, cblabel, True, "ttemp.png")
        input_im = mimage.imread("ttemp.png")
        a.imshow(input_im)
        a.text(0.01,0.92,label,transform=a.transAxes,bbox=props,size=10)
        a.set_xticks([])
        a.set_yticks([])
    os.remove("ttemp.png")

    for ax in f.get_axes():
        for sp in ax.spines.values():
            sp.set_visible(False)

    cb = f.colorbar(raw_im, ax=axes, shrink=0.75)
    cb.ax.tick_params(labelsize=8)
    if cblabel is not None : cb.ax.set_title(cblabel, size=8)

    f.savefig(filename,format="png",dpi=300)

if __name__ == '__main__' :

    parser = argparse.ArgumentParser(description="Plotting 2D maps")
    parser.add_argument('-f', '--file',nargs="+",help="the input file")
    parser.add_argument('--cutoff', type=float, default=0.0, help="the cutoff for filtering")
    parser.add_argument('--max', type=float, help="the maximum value")
    parser.add_argument('--xlabel',help="the x label")
    parser.add_argument('--ylabel',help="the y label")
    parser.add_argument('--cblabel',help="the label for the colorbar")
    parser.add_argument('--extent',type=float,nargs=4,help="the extent (left, right, bottom, top) of the image")
    parser.add_argument('--nobar',action="store_true",help="don't draw a colorbar", default=False)
    parser.add_argument('--multiplot',action="store_true",help="draw more maps in subplots", default=False)
    parser.add_argument('-o','--out',help="the output filename of the multiplot")
    args = parser.parse_args()

    all_mat = [parsing.parse2ndarray(filename) for filename in args.file]

    if not args.multiplot :
        for filename, mat in zip(args.file, all_mat) :
            outname = os.path.splitext(filename)[0]+".png"
            if args.out is not None and len(args.file) == 1 :
                outname = args.out
            _draw_2dmap(mat, args.cutoff, args.max, args.ylabel, args.xlabel, args.extent,
                args.cblabel, args.nobar, outname)
    else :
        _draw_multi_2dmap(all_mat, args.cutoff, args.max, args.ylabel, args.xlabel, args.extent,
            args.cblabel, args.out)
