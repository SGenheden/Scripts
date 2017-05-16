# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to make a Lammps datafile from a PDB file
"""

import sys
import math
import random
import argparse
import os
import copy

import numpy as np

from sgenlib import lammps
from sgenlib import pdb

if __name__ == '__main__' :

  # Command-line input
  parser = argparse.ArgumentParser(description="Make a Lammps data file")
  parser.add_argument('file',help="the PDB or GRO file")
  parser.add_argument('-i','--include',help="the LAMMPS include file")
  parser.add_argument('-o','--out',help="the output prefix")
  parser.add_argument('-b','--box',type=float,nargs=3,help="the box dimensions",default=[0.0,0.0,0.0])
  #parser.add_argument('-a','--atomistic',nargs="+",help="data file(s) for atomistic solutes",default=[])
  parser.add_argument('-c','--converter',help="the dictionary with conversion rules")
  parser.add_argument('-k','--kind',help="the kind of particles", default="cg")
  #parser.add_argument('-p','--pairfunc',help="the pair function for the AA",default="lj/charmm/coul/long")
  args = parser.parse_args()

  # Load a converter
  converter = lammps.Aa2Cg()
  if args.converter is None :
    converter.read(lammps.get_filename("aa2cg.dat")) # The default
  else :
    converter.read(args.converter)

  # Load the force field file
  include = lammps.Includefile(args.include)

  # Create a Datafile and PDBFile
  pdbfile = pdb.PDBFile(args.file) # Input PDB
  takeres = [True for res in pdbfile.residues]
  data = lammps.Datafile()

  # Convert residues
  all_coords = []
  nwat = 0
  moli = 0
  for i,(res,takethis) in enumerate(zip(pdbfile.residues,takeres)) :
    if not takethis : continue
    moli += 1
    res2 = res.resname.strip().lower()
    found = False
    for residue in converter.residues :
      if residue.name == res2 :
        coord = residue.generate_cg(res,moli+1,data,mapping=False)
        all_coords.extend(coord)
        found = True
        break
    # If we could not find a conversion, we will convert the residue to a water bead
    if not found :
      for residue in converter.residues :
        if residue.name == "wat" :
          nwat = nwat + 1
          coord = residue.generate_cg(res,0,data,mapping=False)
          all_coords.extend(coord)

  all_coords = np.array(all_coords)

  print "Minimum of coordinates = %.3f %.3f %.3f"%tuple(all_coords.min(axis=0))
  print "Maximum of coordinates = %.3f %.3f %.3f"%tuple(all_coords.max(axis=0))
  print "Average of coordinates = %.3f %.3f %.3f"%tuple(all_coords.mean(axis=0))

  # Setting correct type for all atoms
  for atom in data.atoms :
    atom.kind = args.kind

  # Settings the correct number of atom and connectivity types
  data.atomtypes = [None]*len(include.masses)
  data.bondtypes = [None]*len(include.bondparams)
  data.angletypes = [None]*len(include.angleparams)
  data.dihedraltypes = [None]*len(include.dihedralparams)

  # Setting the box of the datafile
  if all_coords.mean(axis=0).sum() > 10 : # Checking if center is at origin or not
    data.box = [0.0,0.0,0.0,args.box[0],args.box[1],args.box[2]]
  else :
    data.box = [-args.box[0]/2.0,-args.box[1]/2.0,-args.box[2]/2.0,args.box[0]/2.0,args.box[1]/2.0,args.box[2]/2.0]

  # Write out datafile
  print "Saving LAMMPS data file to %s"%("data."+args.out)
  data.write("data."+args.out)
