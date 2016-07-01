# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to flip the z-coordinate

Examples
--------

"""

import argparse

import numpy as np

from sgenlib import pdb

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Program to rename atoms and residue")
    parser.add_argument('file',help="the PDB file")
    parser.add_argument('-o','--out',help="the output filename")
    args = parser.parse_args()

    pdbin = pdb.PDBFile(args.file)
    pdbin.box[2] *= 2
    for atom in pdbin.atoms:
        new = np.array(atom.xyz, copy=True)
        new[2] = pdbin.box[2]-new[2]
        atom.set_xyz(new)

    pdbin.write(args.out)
