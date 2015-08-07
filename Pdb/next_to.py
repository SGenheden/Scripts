# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to put one PDB-file next to another

If no output file is given the mobile PDB file will be overwritten

Examples
--------
  next_to.py -r pdb1.pdb -m pdb2.pdb
  next_to.py -r pdb1.pdb -m pdb2.pdb -s 0.0 10.0 0.0
  next_to.py -r pdb1.pdb -m pdb2.pdb -s 0.0 10.0 0.0 --fromedge
"""

import argparse
import os
import sys

import numpy as np

from sgenlib import pdb


if __name__ == "__main__":


  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Program to put one PDB-file next to another")
  parser.add_argument('-reference','--reference',help="the reference PDB file")
  parser.add_argument('-mobile','--mobile',help="the mobile PDB file")
  parser.add_argument('-s','--shift',type=float,nargs=3,help="the amount to shift in each direction")
  parser.add_argument('--fromedge',action="store_true",help="if the shift should be from the edge of the mobile structure",default=False)
  parser.add_argument('-o','--out',help="the output filename, default=[overwrite]")
  args = parser.parse_args()

  pdbref = pdb.PDBFile(args.reference)
  pdbmob = pdb.PDBFile(args.mobile)

  # Calculate the delta that will make the centroid of the two structures to overlap
  avref = np.average(pdbref.xyz,axis=0)
  avmob = np.average(pdbmob.xyz,axis=0)
  delta = avref-avmob

  # By command line arguments we can shift the overlay
  if args.shift is None :
    pdbmob.update_xyz(pdbmob.xyz+delta)
  else :
    # If we should place the mobile next to the edge of the reference structure
    if args.fromedge :
      delta2 = pdbref.xyz[:,0].max() - avref
      args.shift = [s+d if abs(s) > 0 else 0 for s,d in zip(args.shift,delta2) ]
    pdbmob.update_xyz(pdbmob.xyz+delta+np.array(args.shift))

  if args.out is None :
    pdbmob.write(args.mobile)
  else :
    pdbmob.write(args.out)

  
