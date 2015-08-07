# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to make input to Gaussian for (R)ESPA calculations

Examples
--------
  pdb2respa_inp.py mol1.pdb mol2.pdb 
  pdb2respa_inp.py mol1.pdb -v ff03 -c -1
  
"""

import sys
import os
import argparse
import time

import numpy as np

from sgenlib import pdb

if __name__ == '__main__' :

  # Command-line input
  parser = argparse.ArgumentParser(description="Making Gaussian input for (R)ESPA calculations")
  parser.add_argument('file',nargs="+",help="the PDB files")
  parser.add_argument('-v','--version',choices=["ff94","ff03"],help="the force field version, can be either ff94 or ff03",default="ff94")
  parser.add_argument('-c','--charge',type=float,help="The net charge of the molecule(s)",default=0)
  parser.add_argument('-p','--processors',type=int,help="The number of processors to use",default=1)
  args = parser.parse_args()

  method = {"ff94" : "HF/6-31G* SCF","ff03" : "B3LYP/cc-pVTZ SCRF"}

  for filename in args.file :
    h,t = os.path.splitext(filename)
    pdbfile = pdb.PDBFile(filename=filename)
    with open("%s_mk.com"%h,"w") as fout : 
      fout.write("%Mem=256MB\n")
      fout.write("%snproc=%d\n"%('%',args.processors))
      fout.write("\n")
      fout.write("#P %s Pop=(Minimal,MK) IOp(6/33=2,6/41=10,6/42=17)\n\n"%method[args.version])
      fout.write("MK ESP on %s, at %s\n"%(filename,time.strftime("%d/%m/%Y")))
      fout.write("\n")
      fout.write("%d 1\n"%args.charge)
      for atom in pdbfile.atoms :
        fout.write("%s %8.3f %8.3f %8.3f\n"%(atom.element(),atom.x,atom.y,atom.z))
      fout.write("\n")
      fout.write("\n")
