# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to make a Gromacs ndx-file from
a configuration file of atom groups and
snapshot
"""

import argparse

from sgenlib import pdb
from sgenlib import groups

if __name__ == "__main__":

  parser = argparse.ArgumentParser(description="Program make density group files")
  parser.add_argument('-f','--file',help="the input gro-file")
  parser.add_argument('-g','--groups',help="the group definitions")
  parser.add_argument('-o','--out',help="the output file",default="density_groups.ndx")
  args = parser.parse_args()

  residue_groups = groups.read_groups(args.groups)

  # Loop over all residues in the input pdb-file
  for residue in pdb.PDBFile(filename=args.file).residues :
    if residue.resname not in residue_groups : continue

    for atom in residue.atoms :
      for group in residue_groups[residue.resname].groups :
        if not group.has_atom(atom.name) : continue
        group.add_serial(atom.name,atom.serial)

  # Check for consistency
  for residue in residue_groups :
    residue_groups[residue].check_consistency()
    print "[%s]"%residue
    for g in residue_groups[residue].groups :
      if g.expressed :
        print "\t%s = %d"%(g.name,len(g.serials[g.atoms[0]]))

  # Write output file
  with open(args.out,"w") as f :
    for residue in residue_groups :
      residue_groups[residue].write_ndx(f)
