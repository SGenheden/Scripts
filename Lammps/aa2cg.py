# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to build a CG or AA/CG dual-resolution
system from an AA system. Converts AA residues
automatically to CG using dictionary of conversions.

Writes out LAMMPS datafiles and inclusion file
for the force field, as well as a PDB-file.
"""

import sys
import math
import random
import argparse
import os
import copy

import numpy as np

import lammps

# Import the PDB module
thispath = os.path.dirname(os.path.abspath(__file__))
oneup = os.path.split(thispath)[0]
sys.path.insert(0,os.path.join(oneup,"Pdb"))
import pdb


def _make_pdbres(coords,atom_names,res_name,pdbfile) :
  """
  Adds a residue + atoms to a PDBFile structure
  """
  res = pdb.Residue()   
  for i,(coord,name) in enumerate(zip(coords,atom_names)) :
    patom = pdb.Atom()
    patom.idx = len(pdbfile.atoms)
    patom.hetatm = False
    patom.serial = len(pdbfile.atoms)+1
    patom.name = name
    patom.residue = len(pdbfile.residues)+1
    patom.resname = res_name
    patom.x = coord[0]
    try :
      patom.y = coord[1]
    except :
      print res_name,len(pdbfile.atoms),coords,coord
    patom.z = coord[2]
    res.append(patom)
    pdbfile.atoms.append(patom)
  pdbfile.residues.append(res)

def _ntypes(array) :
  """
  Find the number of types from a list (could be masses or connectivities)
  """
  n = 0
  for m in array : n = max(n,m.idx)
  return n

def _generate_aa_residue(residue,molidx,resdata,sysdata) :
  """
  Generates an aa residue by copying most of the structure
  from a datafile template, but the coordinates from a PDB residue
  """
  n = len(sysdata.atoms) 
  for i,(ratom,datom) in enumerate(zip(residue.atoms,resdata.atoms)) :
    atom = copy.deepcopy(datom)
    atom.idx = atom.idx + n
    atom.set_xyz(ratom.xyz)
    atom.molecule = molidx
    sysdata.atoms.append(atom)

  for bond in resdata.bonds :  
    b = lammps.Connectivity(record="%d %d %d %d"%(len(sysdata.bonds)+1,bond[0],bond[1]+n,bond[2]+n))
    data.bonds.append(b)

  for angle in resdata.angles :      
    a = lammps.Connectivity(record="%d %d %d %d %d"%(len(sysdata.angles)+1,angle[0],angle[1]+n,angle[2]+n,angle[3]+n))
    data.angles.append(a)

  for dihedral in resdata.dihedrals :      
    a = lammps.Connectivity(record="%d %d %d %d %d %d"%(len(sysdata.dihedrals)+1,dihedral[0],dihedral[1]+n,dihedral[2]+n,dihedral[3]+n,dihedral[4]+n))
    data.dihedrals.append(a)

if __name__ == '__main__' :

  # Command-line input
  parser = argparse.ArgumentParser(description="Converting a AA to CG or AA/CG")
  parser.add_argument('-f','--file',help="the PDB or GRO file")
  parser.add_argument('-i','--include',help="the LAMMPS include file")
  parser.add_argument('-o','--out',help="the output prefix",default="converted")
  parser.add_argument('-b','--box',type=float,nargs=3,help="the box dimensions",default=[0.0,0.0,0.0])
  parser.add_argument('-a','--atomistic',nargs="+",help="data file(s) for atomistic solutes",default=[])
  parser.add_argument('-c','--converter',help="the dictionary with conversion rules")
  args = parser.parse_args()

  # Load a converter
  converter = lammps.Aa2Cg()
  if args.converter is None :
    converter.read(lammps.get_filename("aa2cg.dat")) # The default
  else :
    converter.read(args.converter)

  # Create a Datafile and PDBFile
  pdbfile = pdb.PDBFile(args.file) # Input PDB
  data = lammps.Datafile() 
  pdbout = pdb.PDBFile() # Output PDB

  # Load the force field file
  include = lammps.Includefile(args.include)

  if args.atomistic : 
    natomtypes = _ntypes(include.masses)
    contypes = [_ntypes(include.bondparams), _ntypes(include.angleparams), _ntypes(include.dihedralparams)]
    # Set the functional form of the CG particles if will retain some atomistic molecules
    for pair in include.pair_coeff :
      pair.func = "lj/sf/dipole/sf"
      pair.hybrid = 2

  # Load datafiles for given solutes, will assume these are atomistic
  aa_datafiles = {} 
  for sol in args.atomistic :
    res,filename = sol.split("=")
    res = res.lower()
    aa_datafiles[res] = lammps.Datafile(filename)
    # Extend the force field parameters
    # by extending the inclusion file, we will automatically update the parameter index
    include.extend_from_data(aa_datafiles[res],lj_hybrid=[1,-1],lj_func=["lj/sf/dipole/sf","lj/charmm/coul/long"],ang_func="harmonic")
    # Update the atom and conectivity parameters
    for atom in aa_datafiles[res].atoms : 
      atom.atype = atom.atype + natomtypes
      atom.diameter = 0.0
      atom.density = 1.0
      atom.set_mu([0.0,0.0,0.0])
    conlist = [aa_datafiles[res].bonds,aa_datafiles[res].angles,aa_datafiles[res].dihedrals]
    for cons,ntypes in zip(conlist,contypes) :
      for con in cons : con.param = con.param + ntypes

  # Convert residues
  all_coords = []
  nwat = 0
  for i,res in enumerate(pdbfile.residues) :
    res2 = res.resname.strip().lower()
    found = False
    # If we have an all-atom datafile as a template, keep it as all-atom
    if res2 in aa_datafiles :
      _generate_aa_residue(res,i+1-nwat,aa_datafiles[res2],data)
      coord = res.collect("xyz")
      _make_pdbres(coord,[atom.name for atom in res.atoms],res2,pdbout)
      found = True
    # Otherwise convert it to CG
    else :
      for residue in converter.residues :
        if residue.name == res2 :
          coord = residue.generate_cg(res,i+1,data)
          all_coords.extend(coord)
          _make_pdbres(coord,residue.cg_names,res2,pdbout)
          found = True
          break
    # If we could not find a conversion, we will convert the residue to a water bead
    if not found :
      for residue in converter.residues :
        if residue.name == "wat" :
          nwat = nwat + 1
          coord = residue.generate_cg(res,0,data)
          all_coords.extend(coord)
          _make_pdbres(coord,residue.cg_names,"wat",pdbout)
  all_coords = np.array(all_coords)

  print "Minimum of coordinates = %.3f %.3f %.3f"%tuple(all_coords.min(axis=0))
  print "Maximum of coordinates = %.3f %.3f %.3f"%tuple(all_coords.max(axis=0))
  print "Average of coordinates = %.3f %.3f %.3f"%tuple(all_coords.mean(axis=0))

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

  # Setting correct type for all atoms
  for atom in data.atoms :
    atom.kind = "cg/aa"
 
  # Setting box and adding connectivity to PDB-file
  pdbout.box = args.box
  def add_con(pdbfile,f) :
    for bnd in data.bonds :
      f.write("CONECT%5d%5d\n"%(bnd.atoms[0],bnd.atoms[1]))

  # Write out datafile, pdbfile and force field
  print "Saving LAMMPS data file to %s"%("data."+args.out)
  data.write("data."+args.out)  
  print "Saving PDB file to %s"%(args.out+".pdb")
  pdbout.write(args.out+".pdb",add_extra=add_con)
  if args.atomistic :
    print "Saving LAMMPS inclusion file to %s"%("forcefield."+args.out)
    include.write("forcefield."+args.out)

