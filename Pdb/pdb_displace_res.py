# Author: Samuel Genheden, samuel.genheden@gmail.com

"""
Program to displace a small solute for APR setup
"""

import argparse

from sgenlib import pdb

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description="Script to displace a small solute")
    argparser.add_argument('-f','--file',help="the input structure",default="align_z.pdb")
    argparser.add_argument('-o','--output',help="the output structure file",default="align_z.pdb")
    argparser.add_argument('-r','--resname',help="the residue name",default="MOL")
    argparser.add_argument('-d','--displace',type=float,nargs=+3,help="the displace vector",default=[0.0, 0.0, 0.0])
    args = argparser.parse_args()

    struct = pdb.PDBFile(args.file)
    for atom in struct.atoms :
        if atom.resname.strip() == args.resname :
            atom.set_xyz(atom.xyz+args.displace)
    struct.write(args.output)
