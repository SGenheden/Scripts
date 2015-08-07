# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to modify the atom type of some atoms in a datafile

The --atoms argument is a list of atom specification, which is in format
MOLRANGE:ATOMID

The default output is the input data file with a "_mod" string appended

Examples
  lmp_mod_atype.py data.128dopc_4232wat -a 1-128:9 1-128:15 -t 7 --increase
"""

import argparse

from sgenlib import lammps

if __name__ == '__main__' :

  # Command-line input
  parser = argparse.ArgumentParser(description="Changing atom type of specific atoms")
  parser.add_argument('data',help="the lammps data file")
  parser.add_argument('-o','--out',help="the output file")
  parser.add_argument('-a','--atoms',nargs="+",help="the atoms to change")
  parser.add_argument('-t','--type',type=int,help="the new atom type")
  parser.add_argument('--increase',help="increase number of atom types",action='store_true',default=False)
  parser.add_argument('--decrease',help="decrease number of atom types",action='store_true',default=False)
  args = parser.parse_args()

  # Read datafile and parse into molecules
  datafile = lammps.Datafile(filename=args.data)
  mols = lammps.parse_molecules(datafile)
  
  # Loop over all atom selections and modify atom type 
  for atom in args.atoms :
    molrange,atomi = atom.split(":")
    atomi = int(atomi)
    first,last = map(int,molrange.split("-"))
    print "Changing: ",first,last,atomi
    for m in range(first,last+1) :
      mols[m][atomi-1].atype = args.type

  # Increase or decrease the atomtype list
  if args.increase :
    datafile.atomtypes.append(None)
  elif args.decrease :
    if args.type - 1 < len(datafile.atomtypes) :
      datafile.atomtypes.pop(args.type-1)

  if args.out is None :
    args.out = args.data+"_mod"
  datafile.write(args.out)      
