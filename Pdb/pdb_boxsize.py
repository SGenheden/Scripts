# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to calculate the size of a box that encompass a PDB structure
"""

import argparse

import numpy as np

from sgenlib import pdb

if __name__ == '__main__' :

  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Program to make a box size for a structure")
  parser.add_argument('file',help="the PDB file")
  parser.add_argument('-p','--padding',type=float,nargs=3,help="the padding of each axis",default=[10, 10, 10])
  args = parser.parse_args()

  coord = pdb.PDBFile(filename=args.file).xyz
  av = np.average(coord,axis=0)
  mina = np.min(coord,axis=0)
  maxa = np.max(coord,axis=0)

  lo = mina - args.padding
  hi = maxa + args.padding
  len = hi - lo
  print "Low: %.3f %.3f %.3f"%tuple(lo)
  print "High: %.3f %.3f %.3f"%tuple(hi)
  print "Length: %.3f %.3f %.3f"%tuple(len)
