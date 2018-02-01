# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to perform hydrogen mass repartitioning of a Gromacs topology file

Used in membrane engineering project

Examples
--------
gmx_hmassrepart.py -f dopc.top -o dopc_heavyh.itp
"""
import argparse
import os
import sys

import parmed

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Program to do mass repartitioning")
    parser.add_argument('-f','--file',help="the input top-file")
    parser.add_argument('-o','--out',help="the output top-file",default="hmassrepart.top")
    parser.add_argument('-m','--mass', type=float, help="the new hydrogen mass", default=4.032)
    args = parser.parse_args()

    parm = parmed.load_file(args.file, parametrize=False)
    for i, atom in enumerate(parm.atoms):
        if atom.atomic_number != 1: continue
        heteroatom = None
        heteroidx = 0
        bondeds = list(atom.bond_partners)
        while heteroidx < len(bondeds):
            if bondeds[heteroidx].atomic_number != 1:
                heteroatom = bondeds[heteroidx]
                break
            heteroidx += 1
        transfermass = args.mass - atom.mass
        atom.mass = args.mass
        heteroatom.mass -= transfermass

    if os.path.splitext(args.out)[1] == ".itp" :
        parm.write(args.out, itp=True)
    else :
        parm.write(args.out)
