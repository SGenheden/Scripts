# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to fit an AA structure onto a CG structure by overlaying
the backbone atoms at the BB bead
"""

import argparse
import os
import sys

import numpy as np

from sgenlib import fitting
from sgenlib import pdb

if __name__ == '__main__' :

  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Fitting a AA structure roughly on a CG structure")
  parser.add_argument('-a','--aa',help="all-atom file",default="")
  parser.add_argument('-c','--cg',help="cg file",default="")
  parser.add_argument('-o','--out',help="the output file.",default="out.gro")
  args = parser.parse_args()

  aastruct = pdb.PDBFile()
  aastruct.read(args.aa)

  cgstruct = pdb.PDBFile()
  cgstruct.read(args.cg)

  nres = len(cgstruct.residues)

  # Collect the mobile atoms, i.e. the AA backbone atoms
  mob = np.zeros([nres,3])
  i = 0
  for res in aastruct.residues :
    if res.resname not in pdb.std_aa_names : continue
    cent = np.zeros(3)
    for atom in res.atoms :
      if atom.name.strip() in ["N","CA","C","O"] :
        cent = cent + atom.xyz
    mob[i,:] = cent / 4.0
    i = i + 1

  # Collect the reference atoms, i.e. the BB beads
  ref = np.zeros([nres,3])
  for i,res in enumerate(cgstruct.residues) :
    for atom in res.atoms :
      if atom.name.strip() == "BB" :
        ref[i,:] = atom.xyz
        break

  print "Initial RMSD=",np.sqrt(np.mean(np.sum((mob-ref)**2,axis=1)))

  xyz2 = fitting.dofit(ref,mob,aastruct.xyz)

  for i,c in enumerate(xyz2) :
    aastruct.atoms[i].x = c[0]
    aastruct.atoms[i].y = c[1]
    aastruct.atoms[i].z = c[2]
    aastruct.atoms[i].xyz = c
    
  aastruct.write_gro(args.out)
