# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to replace the coordinates in a datafile with those of a PDB

The default output is the input data file with a "_repxyz" string appended

Example:
  lmp_replace_dataxyz.py -f data.200popc_b2 -p 200popc_b2_pushed.pdb
"""

import os
import sys
import argparse

import numpy as np

from sgenlib import pdb
from sgenlib import lammps

if __name__ == '__main__' :

  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Replace coordinates in LAMMPS datafile with coordinate from a PDB")
  parser.add_argument('-f','--file',help="the name of the datafile ")
  parser.add_argument('-p','--pdb',help="the name of the pdb file ")
  parser.add_argument('-o','--out',help="the name of the output")
  args = parser.parse_args()

  # Read data file and pdb file
  datafile = lammps.Datafile(filename=args.file) 
  pdbfile  = pdb.PDBFile(args.pdb)

  # Replace coordinates
  for da,pa in zip(datafile.atoms,pdbfile.atoms) :
    da.set_xyz(pa.xyz)
 
  # Set box
  datafile.box[:3] = pdbfile.xyz.min(axis=0)-0.1
  datafile.box[3:] = pdbfile.xyz.max(axis=0)+0.1
  
  if args.out is None :
    args.out = args.file+"_repxyz"
  datafile.write(args.out,writeparams=True)
