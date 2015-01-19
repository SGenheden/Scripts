# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to sort the residues in a PDB file
"""

import sys
import argparse
import os

import numpy as np

import pdb

if __name__ == '__main__' :

  # Command-line input
  parser = argparse.ArgumentParser(description="Sort the residue in a PDB file")
  parser.add_argument('-f','--file',help="the PDB file")
  parser.add_argument('-o','--out',help="the output file",default="sorted.pdb")
  parser.add_argument('-r','--residues',nargs="+",help="the residue order",default=[])
  args = parser.parse_args()

  # PDB structures
  pdbin = pdb.PDBFile(filename=args.file)
  pdbout = pdb.PDBFile()

  # Sort the PDB residues in a dictionary
  residues = [r.strip().lower() for r in args.residues]
  resdic = {}
  for residue in pdbin.residues :
    resnam = residue.resname.strip().lower()
    if resnam in residues :
      if resnam not in resdic :
        resdic[resnam] = []  
      resdic[resnam].append(residue)

  # Add the residues sorted in a new PDB object
  for resnam in residues :
    pdbout.extend_residues(resdic[resnam],makecopy=False,dochains=False)
  pdbout.renumber()
  pdbout.write(args.out)
