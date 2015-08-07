# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to calculate the extent of PDB structure
"""

import sys

import numpy as np

from sgenlib import pdb

if __name__ == '__main__' :

  coord = pdb.PDBFile(filename=sys.argv[1]).xyz
  av = np.average(coord,axis=0)
  mina = np.min(coord,axis=0)
  maxa = np.max(coord,axis=0)
  print "Mean: %.3f %.3f %.3f"%(av[0], av[1],av[2])
  print "Min:  %.3f %.3f %.3f"%(mina[0], mina[1],mina[2])
  print "Max:  %.3f %.3f %.3f"%(maxa[0], maxa[1],maxa[2])
  print "Len:  %.3f %.3f %.3f"%(maxa[0]-mina[0], maxa[1]-mina[1],maxa[2]-mina[2])
