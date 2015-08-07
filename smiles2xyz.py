# Author: Samuel Genheden, samuel.genheden@gmail.com

"""
Program to convert a SMILES string to 3D coordinates using
one of two web services. 

Examples
--------
  smiles2xyz.py COH
  smiles2xyz.py -o ethanol.cyx CCOH
"""

import argparse
 
from sgenlib import smiles
   
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Program to convert a SMILES string to 3D-coordinates in xyz format",)
    parser.add_argument('smiles',help="the SMILES string")
    parser.add_argument('-o','--out',help="the output xyz-file, default='mol.xyz'",default="mol.xyz")
    parser.add_argument('-w','--web',choices=["cactus","indiana","both"],help="the name of the web service used to convert the SMILES, should be either 'cactus','indiana', or 'both'",default="both")
    args = parser.parse_args()

    smiles.convert2xyz(args.smiles,args.out,args.web)
