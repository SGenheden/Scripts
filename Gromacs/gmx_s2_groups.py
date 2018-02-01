# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to make a Gromacs ndx-file
for deuterium order parameter calculations
"""


import sys
import os
import argparse

from sgenlib import pdb
from sgenlib import groups

if __name__ == "__main__":


  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Program make S2 group files")
  parser.add_argument('-f','--file',help="the input gro-file")
  parser.add_argument('-g','--groups',help="the group definitions")
  parser.add_argument('-o','--out',help="the output prefix",default="")
  args = parser.parse_args()

  residue_groups = groups.read_groups(args.groups)

  # Loop over all residues in the input pdb-file
  for residue in pdb.PDBFile(filename=args.file).residues :
    if residue.resname not in residue_groups : continue
    for atom in residue.atoms :
      for group in residue_groups[residue.resname].groups :
        if not group.has_atom(atom.name.strip()) : continue
        group.add_serial(atom.name.strip(),atom.serial)

  # Check for consistency
  for residue in residue_groups :
    residue_groups[residue].check_consistency()
    print "[%s]"%residue
    for g in residue_groups[residue].groups :
      if g.expressed :
        print "\t%s = %d"%(g.name,len(g.serials[g.atoms[0]]))

  # Write output file
  for residue in residue_groups :
    for i,g in enumerate(residue_groups[residue].groups,1) :
      if not g.expressed : continue
      with open("%s_chain%d%s.ndx"%(residue.lower(),i,args.out),"w") as f :
        g.write_ndx(f,cutoff=10,split=True)
