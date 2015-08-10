# Author: Samuel Genheden, samuel.genheden@gmail.com

"""
Program to plot 2D maps
"""

import argparse

import numpy as np
import matplotlib.pylab as plt

def _draw_2dmap(mat,filename) :

  f = plt.figure()
  im = f.gca().imshow(mat,cmap=plt.cm.YlGn,origin="lower")
  sel = mat>0
  excluded_im = np.ones([mat.shape[0],mat.shape[1],4])
  excluded_im[sel,3] = 0.0
  f.gca().imshow(excluded_im,origin="lower")
  f.colorbar(im)
  f.savefig(filename,format="png")
  print filename

if __name__ == '__main__' :

  parser = argparse.ArgumentParser(description="Plotting 2D maps")
  parser.add_argument('file',help="the npz file.")
  parser.add_argument('-o','--out',help="a output prefix")
  parser.add_argument('-m','--maps',nargs="+",help="the maps to plot")
  args = parser.parse_args()

  npzfile = np.load(args.file)
  for m in args.maps :
    if m in npzfile.files :
      _draw_2dmap(npzfile[m],args.out+"_"+m+".png")
    else :
      print "Skipping %s"%m
