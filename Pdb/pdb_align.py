# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program align a PDB structure to a specific axis

Examples
--------
  pdb_align.py prot.pdb -a z
"""

import argparse
import os
import sys

import numpy as np

from sgenlib import fitting
from sgenlib import geo
from sgenlib import pdb

if __name__ == "__main__":


  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Program to make align a PDB file to a specific axis")
  parser.add_argument('file',help="the PDB file")
  parser.add_argument('-a','--axis',choices=["x","y","z"],help="the axis to align to",default="z")
  parser.add_argument('-o','--out',help="the output filename")
  args = parser.parse_args()

  pdbfile = pdb.PDBFile(filename=args.file)
  normv = np.zeros(3)
  normv[["x","y","z"].index(args.axis)] = 1.0

  # Collect masses of all residues
  masses = np.zeros(pdbfile.xyz.shape[0])
  iatom = 0
  for residue in pdbfile.residues :
    masses[iatom:iatom+len(residue.atoms)] = residue.collect("masses")
    iatom = iatom+len(residue.atoms)

  # Move protein to center of mass
  center = pdbfile.xyz.mean(axis=0)
  xyz = pdbfile.xyz - center

  # Calculate principl axis and align with user selected axis
  moi = geo.moment_of_inertia(xyz,masses)
  princip = geo.principal_axes(moi)
  rotvec = geo.rotaxis(princip[0,:],normv)
  alpha = geo.angle(princip[0,:],normv)
  rotmat = geo.rotation_matrix(alpha,rotvec)
  xyz = fitting.rotate(xyz,rotmat)+center
  pdbfile.update_xyz(xyz)

  if args.out is None :
    args.out = os.path.splitext(args.file)[0]+"_aligned.pdb"
  pdbfile.write(args.out)
