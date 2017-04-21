# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to rename atoms and residue

Examples
--------
  pdb_rename.py mystruct.pdb -atoms POPC=names_popc.txt TIP3P=names_wat.txt -from Charmm -to Berger -residues TIP3P=SOL
"""

import argparse
import os
import sys

import numpy as np

from sgenlib import pdb

def _readconv(filename,fro,to) :

  lines = open(filename,"r").readlines()
  names = lines[0].strip().split()
  fromidx = names.index(fro)
  toidx = names.index(to)
  conv = {}
  for line in lines[1:] :
    cols = line.strip().split()
    conv[cols[fromidx]] = cols[toidx]
  return conv

if __name__ == "__main__":

  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Program to rename atoms and residue")
  parser.add_argument('file',help="the PDB file")
  parser.add_argument('-o','--out',help="the output filename")
  parser.add_argument('-fro',help="the names in the input")
  parser.add_argument('-to',help="the names in the output")
  parser.add_argument('-atoms',nargs="+",help="the conversion files for atoms")
  parser.add_argument('-residues',nargs="+",help="the conversion for residues")
  args = parser.parse_args()

  pdbin = pdb.PDBFile(args.file)

  if args.atoms is not None :
    convs = {}
    for s in args.atoms :
      res,filename = s.split("=")
      convs[res] = _readconv(filename,args.fro,args.to)
    for residue in pdbin.residues :
      if residue.resname.strip() in convs :
        residue.rename(convs[residue.resname.strip()])

  if args.residues :
    convs = {}
    for s in args.residues :
      f,t = s.split("=")
      convs[f] = t
    for residue in pdbin.residues :
      if residue.resname.strip() in convs :
        new = convs[residue.resname.strip()]
        residue.resname = new
        for atom in residue.atoms : atom.resname = new

  if args.out is None :
    s = os.path.splitext(args.file)
    args.out = s[0]+"_"+args.to.lower()+s[1]
  pdbin.write_gro(args.out)
