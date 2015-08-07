# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to make write out the atom types of heavy atoms

Example:
  lmp_heavy.py forcefield.128dopc_kalp23
"""

import argparse

from sgenlib import lammps

if __name__ == '__main__' :

  # Command-line input
  parser = argparse.ArgumentParser(description="Write out atom types of heavy atoms")
  parser.add_argument('file',help="the lammps include file")
  parser.add_argument('-x','--exclude',type=float,nargs="+",help="exclude these",default=[40.0,90.0,62.0,42.0])
  parser.add_argument('-m','--mass',type=float,help="the mass of hydrogen",default=1.00800)
  args = parser.parse_args()
  
  atypes = []
  for m in lammps.Includefile(args.file).masses :
    if not (m.mass in args.exclude or m.mass == args.mass) :
      atypes.append(m.idx)
  print " ".join(["%d"%i for i in atypes])
