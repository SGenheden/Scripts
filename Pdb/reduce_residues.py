# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program reduce the number of a specific residue, e.g. waters,
by replacing them with a residue at the group centroid.
"""

import argparse
import os

import numpy as np

import pdb

if __name__ == "__main__":

  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Program to reduce residues in a PDB file")
  parser.add_argument('file',help="the PDB file")
  parser.add_argument('-r','--residue',help="the name of the residue to reduce")
  parser.add_argument('-s','--size',type=int,help="the size of the group to reduce")
  parser.add_argument('-o','--out',help="the output filename")
  args = parser.parse_args()

  pdbfile = pdb.PDBFile(filename=args.file)

  taken = [res.resname.strip() != args.residue for res in pdbfile.residues]
  norig = len(pdbfile.residues) - sum(taken)

  if norig % args.size != 0 :
    print "Sorry, but the total number of residues (%d) is not divisble by the group size (%d)."%(norig,args.size)
    print "Cannot proceed. Please change the group size."
    quit()
  else :
    print "Will reduce %d residues into %d."%(norig,norig / args.size)

  ntaken = 0
  i = 0
  while ntaken < norig :
    # Find the new first residue
    first_res = None
    for r,t in zip(pdbfile.residues,taken) :
      if t : continue
      first_res = r
      break
    taken[first_res.idx] = True
    coord = first_res.collect("xyz")
    # Now find args.size - 1 other residues closes to this residue
    dist = []
    res  = []
    for r,t in zip(pdbfile.residues,taken) :
      if t : continue
      dist.append(first_res.distance2(r,pdbfile.xyz))
      res.append(r)
    sortlst = np.asarray(dist).argsort()
    # Hide them and mark them as added
    for isort in sortlst[:args.size-1] :
      res[isort].set_hidden(True)
      taken[res[isort].idx] = True
      coord = coord + res[isort].collect("xyz")
    # Move the first residue to the centroid of the residues
    coord = coord / args.size
    first_res.update_xyz(coord)
    ntaken += args.size

  # Write output
  if args.out is None :
    args.out = os.path.splitext(args.file)[0]+"_reduced.gro"
  pdbfile.write_gro(args.out)
