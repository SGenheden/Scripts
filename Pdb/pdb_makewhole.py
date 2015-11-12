# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to make a PDB structure whole after wrapping in central periodic box

Default name is input with "_whole" appended to the end

Examples
--------
  pdb_makewhole.py prot.pdb
"""

import argparse
import os
import sys

from sgenlib import pbc
from sgenlib import pdb

if __name__ == '__main__' :

  # Command-line input
  parser = argparse.ArgumentParser(description="Making molecules whole in a PDB file")
  parser.add_argument('file',help="the PDB file")
  parser.add_argument('-g','--group',help="group these residue into one entity")
  parser.add_argument('-o','--out',help="the output file")
  args = parser.parse_args()

  pdbfile = pdb.PDBFile(filename=args.file)

  find_group = args.group != None
  if find_group :
    group_first,group_last = args.group.split("-")
    group_first = int(group_first)
    group_last = int(group_last)

  i = 0
  group2 = []
  while i < len(pdbfile.residues) :
    if find_group and group_first == pdbfile.residues[i].serial :
      atoms = []
      while i < len(pdbfile.residues) and pdbfile.residues[i].serial <= group_last :
        atoms.extend(pdbfile.residues[i].atoms)
        i=i+1
      pbc.make_whole(atoms,pdbfile.box)
    else :
      if  len(pdbfile.residues[i].atoms) > 1 :
        pbc.make_whole(pdbfile.residues[i].atoms,pdbfile.box)
      group2.append(pdbfile.residues[i])
      i = i + 1

  if args.out is None :
    args.out = os.path.splitext(args.file)[0]+"_whole.pdb"
  pdbfile.write(args.out,ter=True)
