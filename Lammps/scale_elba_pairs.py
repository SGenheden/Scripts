# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to scale CG-AA pair interactions in LAMMPS inclusion file
"""

import sys
import numpy as np
import argparse

import lammps

if __name__ == '__main__' :

  # Command-line input

  parser = argparse.ArgumentParser(description="Scaling CG-AA interaction in LAMMPS inclusion file")
  parser.add_argument('-f','--file',help="the inclusion file.",default="")
  parser.add_argument('-o','--out',help="the output file.",default="")
  parser.add_argument('-b','--beads',nargs="+",help="the CG atom type groups",default=[])
  parser.add_argument('-s','--scale',type=float,nargs="+",help="the scaling factors",default=[])
  parser.add_argument('-p','--param',nargs="+",help="The parameter to scale, either e(psilon), s(igma), or c(oulomb)",default=[])
  args = parser.parse_args()

  if args.file == "" or args.out == "" :
    print "No input or output files specified. Exiting!"
    quit()

  infile = lammps.Includefile(filename=args.file)

  # Find all mixed pairs
  mixed_pairs = []
  for pair in infile.pair_coeff :
    if pair.comment.find("AA-CG mixed") > -1 or pair.comment.find("Mixed using Lorentz-Berthelot rules") > -1 : mixed_pairs.append(pair)
  lj_same = infile.lj_same_dict()

  # Loop over all the bead specification and change the given parameter
  for beadspec,fac,param in zip(args.beads,args.scale,args.param) :
    beads = [int(b) for b in beadspec.strip().split(",")]
    for pair in mixed_pairs :
      if pair.iatom in beads : 
        if param.lower()[0] == "e" :
          new = fac*np.sqrt(lj_same[pair.iatom].epsilon*lj_same[pair.jatom].epsilon)
          pair.epsilon = new
        elif param.lower()[0] == "s" :
          new = fac*0.5*(lj_same[pair.iatom].sigma+lj_same[pair.jatom].sigma)
          pair.sigma = new
        elif param.lower()[0] == "c" :
          pair.scale = fac
  infile.write(args.out)

 
  
