# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to put one PDB-file next to another

Examples
--------
  next_to.py -r pdb1.pdb -m pdb2.pdb
  next_to.py -r pdb1.pdb -m pdb2.pdb -s 0.0 10.0 0.0
"""

import argparse
import os
import sys

import numpy as np

import pdb


if __name__ == "__main__":


  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Program to put one PDB-file next to another")
  parser.add_argument('-reference','--reference',help="the reference PDB file")
  parser.add_argument('-mobile','--mobile',help="the mobile PDB file")
  parser.add_argument('-s','--shift',type=float,nargs=3,help="the amount to shift in each direction")
  parser.add_argument('-o','--out',help="the output filename, default=[overwrite]")
  args = parser.parse_args()

  pdbref = pdb.PDBFile(args.reference)
  pdbmob = pdb.PDBFile(args.mobile)

  avref = np.average(pdbref.xyz,axis=0)
  avmob = np.average(pdbmob.xyz,axis=0)
  delta = avref-avmob

  if args.shift is None :
    pdbmob.update_xyz(pdbmob.xyz+delta)
  else :
    pdbmob.update_xyz(pdbmob.xyz+delta+np.array(args.shift))


  if args.out is None :
    pdbmob.write(args.mobile)
  else :
    pdbmob.write(args.out)

  
