# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to convert a LAMMPS datafile into a PDB file
using a conversion dictionary or a template PDB file

Example:
  lmp_data2pdb.py data.128dopc_kalp23 -o 128dopc_kalp23_initial.pdb -d 0:wat 1-128:dop 129:kalp23_leap.pdb
"""

import argparse
import os
import sys

import numpy as np

from sgenlib import lammps
from sgenlib import pdb
from sgenlib import pbc

class AtomSelection :
  """
  Class to store atom names, residue names and connectivity for
  a selection of datafile atoms

  Attributes
  ----------
  first : int
    the first molecule in the selection
  last : int
    the last molecule in the selection
  resname : string
    the residue name
  atom_names : list of string
    the atom names
  pdbfile : PDBFile object
    a pdb template if more than one residue
  bonds : list of int tuples
    the connectivity
  """
  def __init__(self,d,converter) :
    ran,nam = d.split(":")

    # Parse the range
    if ran.find("-") > -1 :
      first,last = ran.split("-")
    else :
      first = ran
      last = ran
    self.first = int(first)
    self.last = int(last)

    # Parse the name
    internal_res = {res.name.strip().lower() : res.cg_names for res in converter.residues}
    internal_bonds = {res.name.strip().lower()+"_bonds" : res.bonds for res in converter.residues}
    if nam.strip().lower() in internal_res :
      self.resname = nam
      self.atom_names = internal_res[nam.strip().lower()]
      self.bonds = [(b[1],b[2]) for b in internal_bonds[nam.strip().lower()+"_bonds"]]
    elif os.path.isfile(nam) :
      pdbfile = pdb.PDBFile(filename=nam)
      if len(pdbfile.residues) > 1 :
        self.resname = "FILE!"
        self.atom_names = []
        self.bonds = []
        self.pdbfile = pdbfile
      else :
        self.resname = pdbfile.residues[0].resname
        self.atom_names = [atom.name for atom in pdbfile.residues[0].atoms]
        self.bonds = []
    else :
      raise Exception("Could not find %s in the converter dictionary and it is not a file"%nam)
  def checkout(self,atom) :
    """
    Check if a datafile atom should be parsed by this object
    """
    ok = False
    if self.first == self.last :
      ok = atom.molecule == self.first
    else :
      ok = self.first <= atom.molecule and self.last >= atom.molecule
    return ok


if __name__ == '__main__' :

  # Command-line input
  parser = argparse.ArgumentParser(description="Converting between a lammps data file and a PDB-file")
  parser.add_argument('file',help="the lammps data file")
  parser.add_argument('-o','--out',help="the pdb file")
  parser.add_argument('-d','--dict',help="a dictionary to determine how to parse the data file to a PDB-file",nargs="+",default=["0:W","d:dop"])
  parser.add_argument('-w','--whole',action='store_true',help="turns on making molecules whole",default=False)
  parser.add_argument('-r','--renumber',action='store_true',help="renumber residues",default=False)
  parser.add_argument('-x','--fixed',action='store_true',help="fix residue at the centre of box",default=False)
  parser.add_argument('-b','--boxadjust',action='store_true',help="adjust box dimensions",default=False)
  parser.add_argument('-c','--converter',help="the dictionary with conversion rules")
  args = parser.parse_args()

  # Load a converter
  converter = lammps.Aa2Cg()
  if args.converter is None :
    converter.read(lammps.get_filename("aa2cg.dat")) # The default
  else :
    converter.read(args.converter)

  # Create input and output objects
  data = lammps.Datafile(filename=args.file)
  pdbfile = pdb.PDBFile()

  # Create AtomSelection objects
  selections = [AtomSelection(d,converter) for d in args.dict]

  # Split atoms into molecule
  mols,mollist = lammps.parse_molecules(data,order=True)

  # Loop over all molecules and build the PDB structure
  bonds = []
  group_fix = []
  group_mob = []
  for m in mollist :
    atoms = mols[m]
    selection = None
    for s in selections :
      if s.checkout(atoms[0]) :
        selection = s
        break
    if selection is None :
      raise Exception("Could not find a specified parser for molecule with id %d"%m)

    if selection.resname.lower() == "wat" : # Create a residue for each water molecule
      for a in atoms :
        pdb.make_pdbres([a.xyz],selection.atom_names,selection.resname,pdbfile)
        group_mob.append(pdbfile.residues[-1])
    elif selection.resname == "FILE!" : # Copy an entire PDB object
      selection.pdbfile.update_xyz([a.xyz for a in atoms])
      if args.renumber :
        pdbfile.extend_residues(selection.pdbfile.residues,dochains=False)
      else :
        pdbfile.extend_residues(selection.pdbfile.residues,dochains=False,resnumber=[len(pdbfile.residues)+1])
      n = len(selection.pdbfile.residues)
      for r in pdbfile.residues[-n:] : group_fix.extend(r.atoms)
    else :
      coords = [a.xyz for a in atoms]
      n = len(pdbfile.atoms)
      for b in selection.bonds :
        bonds.append([b[0]+n,b[1]+n])
      pdb.make_pdbres(np.array(coords),selection.atom_names,selection.resname,pdbfile)
      group_mob.append(pdbfile.residues[-1])
  pdbfile.renumber(doresidues=args.renumber)

  # Add box
  pdbfile.box = np.zeros(3)
  pdbfile.box[0] = data.box[3]-data.box[0]
  pdbfile.box[1] = data.box[4]-data.box[1]
  pdbfile.box[2] = data.box[5]-data.box[2]

  # Make each residue whole
  do_fixedgroup = args.fixed and len(group_fix) > 0
  if args.whole or do_fixedgroup :
    print "Making residues whole"
    if do_fixedgroup :
      pbc.make_whole(group_fix,pdbfile.box)
    for residue in pdbfile.residues :
      if len(residue.atoms) == 1 : continue
      pbc.make_whole(residue.atoms,pdbfile.box)
    if do_fixedgroup :
      print "Center fixed group"
      pbc.center(group_fix,group_mob,pdbfile.box)

  if not args.whole and do_fixedgroup :
    print "Putting molecules back in central box"
    pbc.wrap(pdbfile.atoms,pdbfile.box)

  if args.boxadjust :
    print "Adjusting box size if necessary"
    DELTA = 0.1
    xyz = np.array([a.xyz for a in pdbfile.atoms])
    xyz_len = xyz.max(axis=0) - xyz.min(axis=0)
    delta_len = np.abs(xyz_len - pdbfile.box)
    pdbfile.box[delta_len>DELTA] = xyz_len[delta_len>DELTA]+DELTA
    midbox = pdbfile.box / 2.0
    xyz_com = xyz.mean(axis=0)
    delta_com = xyz_com - pdbfile.box / 2.0
    for atom in pdbfile.atoms :
      atom.set_xyz(atom.xyz-delta_com)

  # Write out the PDB file
  if args.out is None :
    s1,s2 = os.path.splitext(args.file)
    args.out = s2[1:]+".pdb"
    print "Writing PDB file: "+args.out

  def add_con(pdbfile,f) :
    for bnd in bonds :
      f.write("CONECT%5d%5d\n"%(bnd[0],bnd[1]))
  pdbfile.write(args.out,ter=True,add_extra=add_con)
