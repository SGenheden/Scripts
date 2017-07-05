# Author: Samuel Genheden, samuel.genheden@gmail.com

"""
Program to convert a SMILES string to 3D coordinates using
one of two web services.

The format of the list when using the --inlist option should be
name1 SMILES1
name2 SMILES2
...

Examples
--------
  smiles2xyz.py COH
  smiles2xyz.py -o ethanol.xyz CCOH
  smiles2xyz.py --inlist smiles.txt
"""

import argparse

from sgenlib import smiles

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Program to convert a SMILES string to 3D-coordinates in xyz format",)
    parser.add_argument('smiles',help="the SMILES string")
    parser.add_argument('-o','--out',help="the output xyz-file, default='mol.xyz'",default="mol.xyz")
    parser.add_argument('-w','--web',choices=["cactus","indiana","both"],help="the name of the web service used to convert the SMILES, should be either 'cactus','indiana', or 'both'",default="both")
    parser.add_argument('--inlist',action="store_true", default=False, help="read names and SMILES from a list")
    args = parser.parse_args()

    if not args.inlist :
        smiles.convert2xyz(args.smiles, args.out, args.web)
    else :
        solutes = [s.strip().split() for s in open(args.smiles,'r').readlines()]
        for solute in solutes :
            smiles.convert2xyz(solute[1], solute[0]+".xyz", args.web)
