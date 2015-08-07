# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to sort the residues in a PDB file

Default name is input with "_sorted" appended to the end

Examples
--------
  sort_residues.py prot.pdb -r SOL WAT CL
"""

import sys
import argparse
import os

import numpy as np

from sgenlib import pdb

if __name__ == '__main__' :

  # Command-line input
  parser = argparse.ArgumentParser(description="Sort the residue in a PDB file")
  parser.add_argument('file',help="the PDB file")
  parser.add_argument('-o','--out',help="the output file")
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
  
  if args.out is None :
    args.out = os.path.splitext(args.file)[0]+"_sorted.pdb"
  pdbout.write(args.out)
