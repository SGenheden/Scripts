# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program plot time series of molecular densities

Examples
--------
gpcr_plot_densityseries.py -f r{1..5}_densities1.npz -o densities_series -d chol -m b2
"""

import os
import argparse

import matplotlib
if not "DISPLAY" in os.environ or os.environ["DISPLAY"] == "" :
  matplotlib.use('Agg')
import matplotlib.pyplot as plt

import gpcr_lib

if __name__ == '__main__' :

  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Plotting pairs of GPCR densities")
  parser.add_argument('-f','--files',nargs="+",help="a list of input files.",default=[])
  parser.add_argument('-o','--out',help="the output png-file.",default="densities_series")
  parser.add_argument('-d','--density',help="the name of the densities.")
  parser.add_argument('-m','--mol',choices=["b2","a2a","b2_a","a2a_a"],help="the protein molecules")
  parser.add_argument('-p','--plot',choices=["dg","average"],help="the density to plot",default="dg")
  parser.add_argument('-r','--replacement',nargs=4,help="the replacement strings",default=["densities2","densities3","densities4","densities5"])
  parser.add_argument('-l','--label',nargs=4,help="the labels of the densities",default=[r'$10 \mu s$',r'$20 \mu s$',r'$30 \mu s$',r'$40 \mu s$'])
  parser.add_argument('--max',type=float)
  parser.add_argument('--min',type=float)
  args = parser.parse_args()

  # Read Xray structure and densities
  xray      = gpcr_lib.load_xray(args.mol)
  densities = [gpcr_lib.standard_read_and_process(args.files,args.density)]
  for r in args.replacement[1:] :
    filenames = [f.replace(args.replacement[0],r) for f in args.files]
    densities.append(gpcr_lib.standard_read_and_process(filenames,args.density))

  if args.plot == "average" :
    for d in densities : d.cutoff_av(cutoff=0.15)

  if args.max is None :
    maxval = max(*[d.max(args.plot) for d in densities])
  else :
    maxval = args.max
  if args.min is None :
    minval = min(*[d.min(args.plot) for d in densities])
  else :
    minval = args.min
  print "Max = %.2f and min = %E"%(maxval,minval)

  # Plot the series
  for i,side in enumerate(["low","upp"],1) : # Loop over sides
    f = plt.figure(i)
    for i,(density,label,number) in enumerate(zip(densities,args.label,["A)","B)","C)","D)"]),1) :
      a = f.add_subplot(2,2,i)
      im = gpcr_lib.plot_density_xray(a,density[side],args.plot,minval,maxval,xray,side,label,number,plotn=False)
    gpcr_lib.draw_colormap(f,im)
    f.savefig("%s_%s_%s.png"%(args.out,args.density,side),format="png",dpi=300)
