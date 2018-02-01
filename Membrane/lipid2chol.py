# Author: Samuel Genheden

"""
Routines to replace MARTINI lipids with MARTINI cholesterol

This module defines a single public function:
replace_lipid

Can be executed from the command line as a stand-alone program
"""

import os

import numpy as np
import numpy.random as random

import pdb
import fitting

lipids = ["POPC","DOPC"]
mapping = {}
mapping["POPC"] = [["C2","D3B"],["C1","C2B"],["R5","C1B"],["R2","C1A"],["R1","GL2"],["ROH","GL1"]]
mapping["DOPC"] = [["C2","D3B"],["C1","C2B"],["R5","C1B"],["R2","C1A"],["R1","GL2"],["ROH","GL1"]]
head_atom = {"POPC" : "PO4", "DOPC" : "PO4"}

def _fit_chol(lipid,cholesterol,indices) :
  """
  Fit a single cholesterol residue onto a lipid residue
  
  Parameters
  ----------
  lipid : Residue
    the lipid molecule reference
  cholesterol : Residue
    the cholesterol that should be fitted
  indices : dictionary of int
    map of atom indices for overlaying

  Returns
  -------
  Residue
    the fitted cholesterol molecule

  Raises
  ------
  ValueError
    if the atom indices is incorrect
  """
  
  # Make a copy of the cholesterol
  cholesterol = pdb.Residue(copy=cholesterol)

  resnam = lipid.resname.upper().strip()
  lipid_indices = indices[resnam]
  chol_indices = indices[resnam+"_CHOL"]
  
  if min(lipid_indices,chol_indices) < 0 or max(lipid_indices) > len(lipid.atoms) or \
     max(chol_indices) > len(cholesterol.atoms) :
       raise ValueError("Bad map indices")
  
  # Extract the mobile coordinates for the fit
  mob = np.zeros([len(chol_indices),3])
  for i,idx in enumerate(chol_indices) :
    mob[i,:] = cholesterol.atoms[idx].xyz

  # Extract the reference coordinates for the fit
  ref = np.zeros([len(lipid_indices),3])
  for i,idx in enumerate(lipid_indices) :
    ref[i,:] = lipid.atoms[idx].xyz  
    
  # Perform the fitting and return the fitted cholesterol  
  xyz = fitting.dofit(ref,mob,cholesterol.collect("xyz"))
  cholesterol.update_xyz(xyz)
  return cholesterol
  

def replace_lipid(membrane,cholesterol,nreplace) :
  """
  Replace lipids with cholesterol
  
  Parameters
  ----------
  membrane : PDBFile
    structure of the membrane, to be modified
  cholesterol : PDBFile
    structure of a template cholesterol molecule
  nreplace : int
    number of lipids to be replace
    
  Returns
  -------
  PDBFile
    the membrane with cholesterol inserted
    
  Raises
  ------
  ValueError
    nreplace must be divisible with 2
  """

  if nreplace % 2 != 0 :
    raise ValueError("Number of replacements must be divisible by 2")
  
  lipids_wat = ["W"]
  lipids_wat.extend(lipids)
  
  # Store away residues
  lipid_res = [r for r in membrane.residues if r.resname.upper().strip() in lipids ]
  wat_res = [r for r in membrane.residues if r.resname.upper().strip() == "W"  ]
  other_res = [r for r in membrane.residues if r.resname.upper().strip() not in lipids_wat ]
  
  # Determine the centre of the bilayer
  zsum = 0
  for res in lipid_res :
    for atom in res.atoms :
      if atom.name.strip().upper() == head_atom[res.resname.upper().strip()] :
        zsum = zsum + atom.z
        break
  zmid = zsum / float(len(lipid_res))
  
  # Determine which lipids are in the lower leaflet
  lower = [False]*len(lipid_res)
  for i,res in enumerate(lipid_res) :
    for atom in res.atoms :
      if atom.name.strip().upper() == head_atom[res.resname.upper().strip()] :
        lower[i] = atom.z < zmid 
        break
  nlower = sum(lower)
  #print "Found a distribution of %d lipids in the lower leaflet and %d lipids in the upper leaflet"%(nlower,len(lipid_res)-nlower)
          
  # Find the indices of the atoms mapping atoms
  indices = {}
  for res in lipid_res :
    resnam = res.resname.upper().strip()
    if resnam in indices : continue
    indices[resnam] = [-1]*len(mapping[resnam])
    for mi,m in enumerate(mapping[resnam]) :
      for i,atom in enumerate(res.atoms) :
        atomnam = atom.name.strip().upper()
        if atomnam == m[1] : 
          indices[resnam][mi] = i
          break
    indices[resnam+"_CHOL"] = [-1]*len(mapping[resnam])
    for mi,m in enumerate(mapping[resnam]) :
      for i,atom in enumerate(cholesterol.residues[0].atoms) :
        atomnam = atom.name.strip().upper()
        if atomnam == m[0] : 
          indices[resnam+"_CHOL"][mi] = i
          break
  
  # Do the random replacement
  chol_res = []
  taken = [False]*len(lipid_res)
  nreplace2 = nreplace / 2
  while len(chol_res) < nreplace2 : # First in the upper leaflet
    probe = np.random.randint(0,len(lipid_res))
    while taken[probe] or lower[probe] : 
      probe = np.random.randint(0,len(lipid_res))
    taken[probe] = True
    chol_res.append(_fit_chol(lipid_res[probe],cholesterol,indices))
  while len(chol_res) < nreplace : # Then in the lower leaflet
    probe = np.random.randint(0,len(lipid_res))
    while taken[probe] or not lower[probe] : 
      probe = np.random.randint(0,len(lipid_res))
    taken[probe] = True
    chol_res.append(_fit_chol(lipid_res[probe],cholesterol,indices))

  # Construct a new PDBFile object and renumber
  new_membrane = pdb.PDBFile()
  new_membrane.extend_residues(other_res,copy=True)
  new_membrane.extend_residues([r for i,r in enumerate(lipid_res) if not taken[i]],copy=True)
  new_membrane.extend_residues(chol_res,copy=False)
  new_membrane.extend_residues(wat_res,copy=True)
  new_membrane.renumber(doatoms=True,doresidues=True)
  new_membrane.box = np.array(membrane.box,copy=True)
  return new_membrane

#
# If this is run from the command-line
#
if __name__ == '__main__' :

  import argparse

  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Program replace a MARTINI lipid with a MARTINI cholesterol")
  parser.add_argument('-f','--file',help="the filename of the input file")
  parser.add_argument('-c','--chol',help="the filename of a cholesterol template file")
  parser.add_argument('-n','--nreplace',nargs="+",type=int,help="the number of replacements to do",default=[2])
  parser.add_argument('--percent',action='store_true',help="turns on input as percentages",default=False)
  args = parser.parse_args()
  
  if args.file is None or args.chol is None :
    print "Nothing to do! Exiting"
    quit()
  
  memstruct = pdb.PDBFile()
  memstruct.read(args.file,gro=True)
  
  choltem = pdb.PDBFile()
  choltem.read(args.chol,gro=True)
  
  head,tail = os.path.splitext(args.file)
  
  if args.percent :
    lipid_res = [r for r in memstruct.residues if r.resname.upper().strip() in lipids ]
    nlipids = len(lipid_res)
  
  for n in args.nreplace :
    nn = n
    if args.percent :
      nn = np.floor(nlipids*n/100.0)
    if nn % 2 != 0 :
      nn = nn - 1
      print "%d is not divisible by two so decreasing the number of replacements with one"%n
    newstruct = replace_lipid(memstruct,choltem,nn)
    filename = "%s_%d_chol.gro"%(head,nn)  
    print "Writing new structure file to: %s\n"%filename
    newstruct.write_gro(filename)
      
