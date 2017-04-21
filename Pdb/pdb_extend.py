# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to extend the box of the structure in PDB or GRO format

Example:
  pdb_extend.py pdbin.pdb -e 5 5 -o pdbout.pdb
"""

import argparse

import numpy as np

from sgenlib import pdb

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Program to extend box")
    parser.add_argument('-f','--file',help="the PDB file")
    parser.add_argument('-o','--out',help="the output filename")
    parser.add_argument('-e','--extend',type=float,nargs=6,help="the size to extend the box up and below")
    args = parser.parse_args()

    pdbfile = pdb.PDBFile(filename=args.file)
    for atom in pdbfile.atoms :
        atom.xyz += np.asarray(args.extend[:3])
        atom.x += args.extend[0]
        atom.y += args.extend[1]
        atom.z += args.extend[2]

    pdbfile.box += np.asarray(args.extend[:3])+np.asarray(args.extend[3:])
    if args.out is None :
        pdbfile.write(args.file)
    else :
        pdbfile.write(args.out)
