# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to do hydrogen mass repartition on a lammps datafile

Will increase the mass of the hydrogen and decrease mass of the heavy atoms bonded to
them. 

The default output is the input data file with a "_heavyh" string appended

Examples
  lmp_massrepart.py data.kalp23 -f 3
"""

import argparse
import copy

from sgenlib import lammps

if __name__ == '__main__' :

  # Command-line input
  parser = argparse.ArgumentParser(description="Doing mass repartitioning on a LAMMPS datafile")
  parser.add_argument('data',help="the lammps data file")
  parser.add_argument('-o','--out',help="the name of the output")
  parser.add_argument('-f','--factor',type=int,help="how much the hydrogen mass should be increased",default=3)
  args = parser.parse_args()
  
  datafile = lammps.Datafile(filename=args.data)
  
  # Find hydrogen atoms in the datafile
  hmass = 1.00800
  hydrogens = []
  totalmass = 0.0
  for atom in datafile.atoms :
    totalmass += datafile.atomtypes[atom.atype-1].mass
    if abs(datafile.atomtypes[atom.atype-1].mass - hmass) < 0.01 :
      hydrogens.append(atom)
  print "Initial mass = %.3f"%totalmass

  # Find atoms bonded to hydrogen atoms and list their type
  bonded_to = []
  hbonds = [0]*len(datafile.atoms)
  for hydrogen in hydrogens :
    for bond in datafile.bonds : # Looping over all bonds in the system
      pairidx = -1
      # See if this is bond to the hydrogen atom
      if hydrogen.idx == bond.atoms[0] : pairidx = bond.atoms[1]
      elif hydrogen.idx == bond.atoms[1] : pairidx = bond.atoms[0]
      # If we find a bond to the hydrogen...
      if pairidx > -1 :
        bonded_to.append(pairidx) # Store away that atom
        hbonds[pairidx-1] += 1 # Increase the number of bonds to this atom
        break # Escape this loop since we assume hydrogen is only bonded to one other atom
  
  # Create a unique list of atom type + hydrogen bonds
  new_types = []
  for i,atom in enumerate(datafile.atoms) :
    if hbonds[i] > 0 : new_types.append((atom.atype,hbonds[i]))
  new_types = list(set(new_types))
  new_types.sort(key=lambda newtype : (newtype[0],newtype[1]))
  print "Will create new types for:"
  for t in new_types : print "Type %2d with %d bonds to hydrogens"%(t[0],t[1])
  
  # Create new atom types and pair coefficients
  noldtype = len(datafile.atomtypes)
  for ti,new_type in enumerate(new_types) :
    atype,nbonds = new_type
    
    # New atom types
    for atomtype in datafile.atomtypes :
      if atomtype.idx == atype :
        datafile.atomtypes.append(copy.deepcopy(atomtype))        
        break    
    datafile.atomtypes[-1].idx = len(datafile.atomtypes)
    datafile.atomtypes[-1].comment = "# New type, mass decreased from %.5f, %d h-bonds"%(datafile.atomtypes[-1].mass,nbonds)
    datafile.atomtypes[-1].mass -= hmass*(args.factor-1)*nbonds

    # New pair types
    if len(datafile.pairtypes) > 0 :
      new_pairs = []
      for pairtype in datafile.pairtypes :
        if atype in [pairtype.iatom,pairtype.jatom] :
          new_pairs.append(copy.deepcopy(pairtype)) 
          new_pairs[-1].comment += " copy of %d,%d"%(new_pairs[-1].iatom,new_pairs[-1].jatom)
          if atype == pairtype.iatom :
            new_pairs[-1].iatom = noldtype+ti+1
          if atype == pairtype.jatom :
            new_pairs[-1].jatom = noldtype+ti+1
          if new_pairs[-1].jatom < new_pairs[-1].iatom :
            temp = new_pairs[-1].jatom
            new_pairs[-1].jatom = new_pairs[-1].iatom
            new_pairs[-1].iatom = temp
          if pairtype.iatom == pairtype.jatom :
            new_pairs.append(copy.deepcopy(pairtype)) 
            new_pairs[-1].comment += " copy of %d,%d"%(new_pairs[-1].iatom,new_pairs[-1].jatom)
            new_pairs[-1].jatom = noldtype+ti+1
      datafile.pairtypes.extend(new_pairs)
  if len(datafile.pairtypes) > 0 :
    datafile.pairtypes.sort(key=lambda pair: (pair.iatom,pair.jatom)) 

  # Increas mass of hydrogen atom types
  for atomtype in datafile.atomtypes :
    if abs(atomtype.mass - hmass) < 0.01 : 
      if args.factor > 0 :
        atomtype.mass = hmass*args.factor
      atomtype.comment = "# Mass increased from %.5f"%hmass
  
  # Change the atom type of heavy atoms
  for i,atom in enumerate(datafile.atoms) :
    if hbonds[i] == 0 : continue
    for j,new_type in enumerate(new_types) :
      if atom.atype == new_type[0] and hbonds[i] == new_type[1] : atom.atype = noldtype + j + 1 

  totalmass = 0.0
  for atom in datafile.atoms :
    totalmass += datafile.atomtypes[atom.atype-1].mass
  print "Final mass = %.3f"%totalmass
  
  if args.out is None :
    args.out = args.data+"_heavyh"
  datafile.write(args.out,writeparams=True)
