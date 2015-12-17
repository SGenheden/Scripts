# Author: Samuel Genheden, samuel.genheden@gmail.com

"""
Program to make 2D scatter plots of data
"""

import argparse

import numpy as np
import matplotlib.pylab as plt

from sgenlib import colors
from sgenlib import parsing
from sgenlib import series

if __name__ == '__main__' :

  parser = argparse.ArgumentParser(description="Plot xy data")
  parser.add_argument('-f','--file',nargs="+",help="the data files")
  parser.add_argument('-l','--label',nargs="+",help="the labels")
  parser.add_argument('-o','--out',help="the output filename",default="xy.png")
  parser.add_argument('-s','--skip',type=float,help="how much to discard",default=0.0)
  parser.add_argument('--xcol',type=int,help="the column of the x-data",default=0)
  parser.add_argument('--ycol',type=int,help="the column of the x-data",default=1)
  parser.add_argument('--xfactor',type=float,default=1.0,help="multiplicator factor for x")
  parser.add_argument('--yfactor',type=float,default=1.0,help="multiplicator factor for y")
  parser.add_argument('--xlabel',help="the x label")
  parser.add_argument('--ylabel',help="the y label")
  parser.add_argument('--xmax',type=float,help="the maximum of the x-axis")
  parser.add_argument('--ymax',type=float,help="the maximum of the y-axis")
  parser.add_argument('--equil',action="store_true",help="estimate equilibration",default=False)
  args = parser.parse_args()

  f = plt.figure(figsize=(3.33,2.5))
  a = f.add_axes((0.2,0.15,0.75,0.65))

  for i,(filename,label) in enumerate(zip(args.file,args.label)) :
    data = parsing.parse2ndarray(filename)
    if args.skip < 1 :
      n = int(np.floor(args.skip * data.shape[0]))
    else :
      n = args.skip
    data = data[n:,:]
    x = data[n:,args.xcol]*args.xfactor
    y = data[n:,args.ycol]*args.yfactor
    a.plot(x,y,colors.style(i),color=colors.color(i),label=label)
    if args.equil :
        equili = series.find_equilibration(x,y)
        print "%s equilibrated at %.3f"%(label,x[equili])

  if len(args.file) > 1 :
      a.legend(loc='best', fancybox=True, framealpha=0.5,fontsize=8,labelspacing=0.20)
  #if len(args.file) > 1 : a.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0.0,fontsize=8)
  for tick in a.xaxis.get_major_ticks() :
    tick.label.set_fontsize(8)
  for tick in a.yaxis.get_major_ticks() :
    tick.label.set_fontsize(8)
  if args.xlabel is not None : a.set_xlabel(args.xlabel,fontsize=8)
  if args.ylabel is not None : a.set_ylabel(args.ylabel,fontsize=8)
  if args.xmax is not None : a.set_xlim([0,args.xmax])
  if args.ymax is not None : a.set_ylim([0,args.ymax])
  f.savefig(args.out,format="png",dpi=300)
