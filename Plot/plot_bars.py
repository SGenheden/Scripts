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
  parser.add_argument('-f','--file',nargs="+",help="the bar plots")
  parser.add_argument('-l','--label',nargs="+",help="the labels")
  parser.add_argument('-o','--out',help="the output filename",default="bar.png")
  parser.add_argument('--xlabel',nargs="+",help="the x label")
  parser.add_argument('--ylabel',help="the y label")
  args = parser.parse_args()

  f = plt.figure()
  a = f.gca()
  for i,(filename,label) in enumerate(zip(args.file,args.label)) :
    data = parsing.parse2ndarray(filename)
    left = np.arange(data.shape[0]) + 0.38 + i*0.75/float(len(args.file))
    a.bar(left,data[:,0],yerr=data[:,1],width=1.0/float(len(args.file))*0.75,color=colors.color(i),label=label,error_kw=dict(ecolor='black'))

  if len(args.file) > 1 :
    a.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=4, mode="expand", borderaxespad=0.0,fontsize=11)
  for tick in a.xaxis.get_major_ticks() :
    tick.label.set_fontsize(8)
  for tick in a.yaxis.get_major_ticks() :
    tick.label.set_fontsize(8)
  if args.xlabel is not None : a.set_xlabel(args.xlabel,fontsize=8)
  if args.ylabel is not None : a.set_ylabel(args.ylabel,fontsize=8)
  f.savefig(args.out,format="png",dpi=300)
