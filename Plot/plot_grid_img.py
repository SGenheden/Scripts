# Author: Samuel Genheden, samuel.genheden@gmail.com

"""
Program to plot figures on a rectangular grid

The size of the individual figures is specified with the -s
argument (width + height) and should be in inches

The configuration of the grid is specified with the -c argument,
number of rows followed by the number of columns

Examples:
  plot_grid_img.py -f a.png b.png -c 2 1 -s 4.0 6.0
  plot_grid_img.py -f a.png b.png c.png d.png -c 4 1 -s 4.0 6.0
  plot_grid_img.py -f a.png b.png c.png d.png -c 2 2 -s 4.0 6.0
                     -l black red green purple --borders
"""

import argparse
import string
from collections import namedtuple

import numpy as np
import matplotlib.pylab as plt
import matplotlib.image as mimage
import matplotlib.gridspec as gridspec

if __name__ == '__main__' :

    parser = argparse.ArgumentParser(description="Plot images on a grid")
    parser.add_argument('-f','--file',nargs="+",help="the images")
    parser.add_argument('-c','--config',type=int,nargs=2,help="the configuration of the file")
    parser.add_argument('-s','--size',type=float,nargs=2,help="the size of the figure", default=[3.33, 2.5])
    parser.add_argument('-l','--labels',nargs="+",help="the labels")
    parser.add_argument('-o','--out',help="the output filename",default="gridimg.png")
    parser.add_argument('--borders',action="store_true",help="plot border between the figures",default=False)
    args = parser.parse_args()

    if args.config is None:
        args.config = [len(args.file),1]

    WidthHeight = namedtuple("with_height",["width","height"])
    config = WidthHeight(width=args.config[1],height=args.config[0])
    size = WidthHeight(width=args.size[0],height=args.size[1])

    if config.height*config.width < len(args.file) :
        print "Configuration is incompatible with the numbers of specified files"
        quit()

    f = plt.figure(figsize=(config.width*size.width,config.height*size.height),dpi=300)
    f.subplots_adjust(wspace=0.0,hspace=0,left=0.00,right=1.0,top=1.0,bottom=0.00)

    labels = [c+")" for c in string.ascii_uppercase[:len(args.file)]]
    if args.labels is not None and len(args.labels) == len(args.file) :
        labels = ["%s %s"%(lbl0,albl) for lbl0,albl in zip(labels,args.labels)]

    props = dict(facecolor='white',edgecolor='white', alpha=1.0)
    for i,(filename,label) in enumerate(zip(args.file,labels),1) :
        ax = f.add_subplot(config.height,config.width,i)
        im = mimage.imread(filename)
        ax.imshow(im)
        ax.text(0.01,0.92,label,transform=ax.transAxes,bbox=props,size=10)
        ax.set_xticks([])
        ax.set_yticks([])

    for ax in f.get_axes():
        for sp in ax.spines.values():
            sp.set_visible(False)
        if args.borders :
            if not ax.is_first_row():
                ax.spines['top'].set_visible(True)
            if not ax.is_last_row():
                ax.spines['bottom'].set_visible(True)
            if not ax.is_first_col():
                ax.spines['left'].set_visible(True)
            if not ax.is_last_col():
                ax.spines['right'].set_visible(True)

    f.savefig(args.out,dpi=300,format="png")
