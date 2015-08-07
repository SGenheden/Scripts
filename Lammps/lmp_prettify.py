# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to make a datafile prettier

If the datafile was created from a restart, the atoms are for instance not in order, 
this fixes that

The default output is the input data file with a "_pretty" string appended

Examples:
  lmp_prettify.py data.128dmpc_kalp23_pushed
  lmp_prettify.py data.128dmpc_kalp23_pushed --writeparams
"""

import argparse

from sgenlib import lammps

if __name__ == '__main__' :

  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Prettify a lammps datafile")
  parser.add_argument('file',help="the name of the datafile ")
  parser.add_argument('-o','--out',help="the name of the output")
  parser.add_argument('-p','--writeparams',action='store_true',help="write force field parameters",default=False)
  args = parser.parse_args()

  data = lammps.Datafile(args.file)
  data.sort()
  for atom in data.atoms :
    atom.ix = None
    atom.iy = None
    atom.iz = None
 
  if not args.writeparams :
    data.atomtypes = [None]*len(data.atomtypes)
    data.bondtypes = [None]*len(data.bondtypes)
    data.angletypes = [None]*len(data.angletypes)
    data.dihedraltypes = [None]*len(data.dihedraltypes)

  if args.out is None :
    args.out = args.file+"_pretty"
  data.write(args.out,writeparams=args.writeparams)
