# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to plot molecular densities

The program takes as input two sets of files and average
the densities in the two sets, creating two averaged densities.

The densities are created by gpcr_mdanal.py

These densities are then plotted next to each other, for easy comparison, along
with a representation of the backbone X-ray structure

Up to five different pictures are produced:
1) The intracellular (low) densities
2) The extracellular (upp) densities
3) Densities in both leaflets
4) Densities in the middle of the bilayer, if such a density is available
5) The difference between the two average densities, if this is specified on the command line

Examples
--------
gpcr_plot_density.py -f1 Inactive/r{1..5}_densities.npz -o density_b2_inactive
                     -d chol2 popc -l chol popc -m b2 --diff
gpcr_plot_density.py -f1 Inactive/r{1..5}_densities.npz -f2 Active/r{1..5}_densities.npz
                     -o density_b2 -d chol2 chol2 -l inactive active -m b2 b2_a

The first compares the cholesterol and POPC densities in inactive B2
The second compares the cholesterol density in inactive to active B2
(chol2 is the density of cholesterol, excluding buried cholesterols)
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
  parser.add_argument('-f1','--files1',nargs='+',help="a list of input files.",default=[])
  parser.add_argument('-f2','--files2',nargs='+',help="a list of input files.")
  parser.add_argument('-o','--out',help="the output png-file.",default="densities")
  parser.add_argument('-d','--density',nargs=2,help="the names of the densities.",default=[])
  parser.add_argument('-m','--mol',nargs="+",choices=["b2","a2a","b2_a","a2a_a"],help="the protein molecules",default=[])
  parser.add_argument('-l','--label',nargs=2,help="labels for the densities.",default=[])
  parser.add_argument('-p','--plot',choices=["dg","average"],help="the density to plot",default="dg")
  parser.add_argument('-c','--countnam',nargs="+",help="the name of the count array")
  parser.add_argument('--diff',help="compute the difference map",action="store_true",default=False)
  parser.add_argument('--max',type=float)
  parser.add_argument('--min',type=float)
  args = parser.parse_args()

  # If we are going to process different densities from the same files,
  # files2 does not have to be specified
  if args.files2 is None : args.files2 = args.files1

  # Make the count name array
  countnam = [d for d in args.density]
  if args.countnam is not None :
    for cn in args.countnam :
      idx,nam = cn.split(":")
      countnam[int(idx)-1] = nam

  # If user specified only one mol, assume the same for the other densities
  if len(args.mol) == 1 : args.mol.append(args.mol[0])

  # Make up pre-defined labels
  labels = [args.label[0]+"/"+gpcr_lib.side_name["low"]]
  labels.append(args.label[1]+"/"+gpcr_lib.side_name["low"])
  labels.append(args.label[0]+"/"+gpcr_lib.side_name["upp"])
  labels.append(args.label[1]+"/"+gpcr_lib.side_name["upp"])
  labels.append(args.label[0]+"/mid("+gpcr_lib.side_name["low"]+")")
  labels.append(args.label[0]+"/mid("+gpcr_lib.side_name["upp"]+")")
  labels.append(args.label[1]+"/mid("+gpcr_lib.side_name["low"]+")")
  labels.append(args.label[1]+"/mid("+gpcr_lib.side_name["upp"]+")")

  # Read the densities and Xray structures from file
  densities1 = gpcr_lib.standard_read_and_process(args.files1,args.density[0],countnam[0])
  xray1      = gpcr_lib.load_xray(args.mol[0])
  densities2 = gpcr_lib.standard_read_and_process(args.files2,args.density[1],countnam[1])
  xray2      = gpcr_lib.load_xray(args.mol[1])

  # Write out the processed densities to file for debugging
  densities1.write(args.plot,args.out+"_savez1.npz")
  densities2.write(args.plot,args.out+"_savez2.npz")

  # Introduce a cut-off when plotting averages
  if args.plot == "average" :
    densities1.cutoff_av(cutoff=0.15)
    densities2.cutoff_av(cutoff=0.15)

  # Find the maximum and minimum values
  if args.max is None :
    maxval = max(densities1.max(args.plot),densities2.max(args.plot))
  else :
    maxval = args.max
  if args.min is None :
    minval = min(densities1.min(args.plot),densities2.min(args.plot))
  else :
    minval = args.min
  print "Max = %.2f and min = %E"%(maxval,minval)

  # Convenient arrays to loop over
  # we will use them below to make tighter code
  densities = [densities1,densities2]
  xray = [xray1,xray2]
  sides = ["low","upp"]

  # Lower leaflets
  fig_low = plt.figure(1)
  for i in range(2) : # Plot each density in its own subplot
    a = fig_low.add_subplot(1,2,i+1)
    im = gpcr_lib.plot_density_xray(a,densities[i].low,args.plot,minval,maxval,xray[i],"low",labels[i])
  gpcr_lib.draw_colormap(fig_low,im)
  fig_low.savefig(args.out+"_low.png",format="png",dpi=300)

  # Upper leaflet
  fig_upp = plt.figure(2)
  for i in range(2) :
    a = fig_upp.add_subplot(1,2,i+1)
    im = gpcr_lib.plot_density_xray(a,densities[i].upp,args.plot,minval,maxval,xray[i],"upp",labels[i+2])
  gpcr_lib.draw_colormap(fig_upp,im)
  fig_upp.savefig(args.out+"_upp.png",format="png",dpi=300)

  # Lower/Upper leaflet
  """fig_lowupp = plt.figure(11)
  a = fig_lowupp.add_subplot(1,2,1)
  im = gpcr_lib.plot_density_xray(a,densities[0].low,args.plot,minval,maxval,xray[i],"low",labels[0],plotn=False)
  a = fig_lowupp.add_subplot(1,2,2)
  im = gpcr_lib.plot_density_xray(a,densities[1].upp,args.plot,minval,maxval,xray[i],"upp",labels[3],plotn=False)
  gpcr_lib.draw_colormap(fig_lowupp,im)
  fig_lowupp.savefig(args.out+"_lowupp.png",format="png",dpi=300)"""

  # Both leaflets
  fig_both = plt.figure(10)
  for i,(dxi,si,label,number) in enumerate(zip([0,1,0,1],[0,0,1,1,],labels[:4],["A)","B)","C)","D)"]),1) :
    # dxi is the index into the densities and xray array,
    #    dxi = 0 gives densities/xray for first group, dxi = 1 gives the second group
    # si is the index into the sides array
    #    si = 0 gives "low", si = 1 gives "upp"
    # these two indices are necessary since the arrays are only 2 items longs, and
    # we don't want to repeat items in the arrays or make different arrays
    # for the mid densities below that uses another order
    a = fig_both.add_subplot(2,2,i)
    im = gpcr_lib.plot_density_xray(a,densities[dxi][sides[si]],args.plot,minval,maxval,xray[dxi],sides[si],label,number=number,plotn=False)
  gpcr_lib.draw_colormap(fig_both,im)
  fig_both.savefig(args.out+"_both.png",format="png",dpi=300)

  # Middle leaflet
  if densities1.mid is not None or densities2.mid is not None :
    fig_mid = plt.figure(3)
    for i,(dxi,si,label,number) in enumerate(zip([0,0,1,1],[0,1,0,1],labels[4:],["A)","B)","C)","D)","E)","F)"]),1) :
      if densities[dxi].mid is not None :
        a = fig_mid.add_subplot(2,2,i)
        im = gpcr_lib.plot_density_xray(a,densities[dxi].mid,args.plot,minval,maxval,xray[dxi],sides[si],"",number="",plotn=False,drawchol=False)
        #im = gpcr_lib.plot_density_xray(a,densities[dxi].mid,args.plot,minval,maxval,xray[dxi],sides[si],label,number=number,plotn=False,drawchol=False)
    gpcr_lib.draw_colormap(fig_mid,im,text=r'$\ln[\rho/2\rho_0]$')
    fig_mid.savefig(args.out+"_mid.png",format="png",dpi=300)

  # Difference
  if args.diff :
    densities1.subtract(args.plot,densities2)
    fig_diff = plt.figure(4)
    for i,(side,label) in enumerate(zip(["low","upp"],[gpcr_lib.side_name["low"],gpcr_lib.side_name["upp"]]),1) :
      a = fig_diff.add_subplot(1,2,i)
      im = gpcr_lib.plot_density_xray(a,densities1[side],args.plot,minval,maxval,xray1,side,label,plotn=False)
    gpcr_lib.draw_colormap(fig_diff,im)
    fig_diff.savefig(args.out+"_diff.png",format="png",dpi=300)
