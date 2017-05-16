# Author: Samuel Genheden samuel.genheden@gmail.com

import argparse

import numpy as np
import numpy.random as random

from sgenlib import pdb

if __name__ == '__main__' :

    parser = argparse.ArgumentParser(description="Split a membrane, making a hole in the middle")
    parser.add_argument('-b','--box',help="the membrane box")
    parser.add_argument('-o','--out',help="the output file",default="splitted.gro")
    parser.add_argument('-p','--patom',help="the name of the phosphate atom", default="P")
    parser.add_argument('-d','--displacement',type=float,help="the displacement length",default="10.0")
    args = parser.parse_args()

    boxfile = pdb.PDBFile(filename=args.box)

    # Calculate the center of membrane
    memz = []
    for atom in boxfile.atoms :
        if atom.name.strip() == args.patom :
            memz.append(atom.z)
    memcent = np.asarray(memz).mean()
    print "Membrane center = %.3f"%memcent

    for residue in boxfile.residues :
        com = residue.collect("centerofmass")
        if com[2] > memcent :
            for atom in residue.atoms :
                new_xyz = atom.xyz + [0,0,args.displacement]
                atom.set_xyz(new_xyz)

    boxfile.box += [0,0,args.displacement]
    boxfile.write(args.out)
