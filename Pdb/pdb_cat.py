# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to concatenate structure in PDB or GRO format

Example:
  pdb_cat.py pdb1.pdb pdb2.pdb pdb3.pdb -o pdb_all.pdb
"""

import argparse
import os

from sgenlib import pdb

if __name__ == "__main__":

    # Setup a parser of the command-line arguments
    parser = argparse.ArgumentParser(description="Program to concatenate structures")
    parser.add_argument('file',nargs="+",help="the PDB file")
    parser.add_argument('-o','--out',help="the output filename")
    args = parser.parse_args()

    pdbfile = pdb.PDBFile(filename=args.file[0], renumber=False)
    for filename in args.file[1:]:
        pdbfile2 = pdb.PDBFile(filename, renumber=False)
        pdbfile.extend_residues(pdbfile2.residues)
    pdbfile.write(args.out)
