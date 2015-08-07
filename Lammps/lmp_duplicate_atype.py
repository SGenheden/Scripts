# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to duplicate atom types in a include file

The default output is the input include file with a "_mod" string appended

Example:
  lmp_duplicate_atype.py forcefield.elba
"""

import argparse

from sgenlib import lammps

if __name__ == '__main__' :

  # Command-line input
  parser = argparse.ArgumentParser(description="Duplicating atom type in Lammps include file")
  parser.add_argument('file',help="the lammps include file")
  parser.add_argument('-o','--out',help="the output file")
  parser.add_argument('-t','--type',type=int,help="the type to duplicate")
  args = parser.parse_args()
  
  ifile = lammps.Includefile(args.file)

  # Duplicate masses and pair coefficients
  newmass,newpairs = ifile.duplicate_type(args.type)
  newmass.comment = "# duplicate of %d"%(args.type)  
  ifile.masses.append(newmass)
  ifile.pair_coeff.extend(newpairs)
  ifile.pair_coeff.sort(cmp=lammps.comp_pair)  

  if args.out is None :
    args.out = args.file+"_mod"
  ifile.write(args.out)
