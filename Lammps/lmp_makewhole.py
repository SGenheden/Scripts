# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to make the molecules in a datafile whole

The default output is the input data file with a "_whole" string appended

Example:
  lmp_makewhole.py data.128dopc_4232wat
"""

import argparse
import os
import sys

import numpy as np

from sgenlib import lammps
from sgenlib import pdb
from sgenlib import pbc

if __name__ == '__main__' :

  # Command-line input
  parser = argparse.ArgumentParser(description="Making molecules whole in a LAMMPS datafile")
  parser.add_argument('file',help="the lammps data file")
  parser.add_argument('-o','--out',help="the output file")
  args = parser.parse_args()

  data = lammps.Datafile(filename=args.file)

  # Create a box with the length of the datafile box
  box = np.zeros(3)
  box[0] = data.box[3]-data.box[0]
  box[1] = data.box[4]-data.box[1]
  box[2] = data.box[5]-data.box[2]

  mols = lammps.parse_molecules(data)
  for atom in data.atoms :
    if atom.ix != None :
      atom.ix = None
      atom.iy = None
      atom.iz = None

  for mol in mols :
    # Assume this is water, have to make this better in the future
    if len(mols[mol]) < 2 or len(mols[mol]) > 200 : continue
    pbc.make_whole(mols[mol],box) 

  # Adjust box
  xyz = np.array([a.xyz for a in data.atoms])
  data.box[:3] = xyz.min(axis=0)
  data.box[3:] = xyz.max(axis=0)

  if args.out is None :
    args.out = args.file+"_whole"
  data.write(args.out)
