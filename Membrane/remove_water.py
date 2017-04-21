# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to reduce the number of water in a membrane, from
e.g. 50 to 40 waters per lipid
"""

import argparse

import numpy as np

from sgenlib import pdb

if __name__ == '__main__' :

    parser = argparse.ArgumentParser(description="Remove water from mebrane")
    parser.add_argument('-f','--file',help="the input file")
    parser.add_argument('-o','--out',help="the output file",default="removed.gro")
    parser.add_argument('-p','--patom',help="the name of the phosphate atom", default="P")
    parser.add_argument('-w','--watom',help="the name of the water atom", default="OH2")
    parser.add_argument('-n','--n',type=int,nargs=2,help="the number of water to remove in each monolayer")
    parser.add_argument('--nobox',action="store_true",help="do not modify the box",default=False)
    args = parser.parse_args()


    struct = pdb.PDBFile(filename=args.file)

    # Calculate the center of membrane
    memz = []
    for atom in struct.atoms :
        if atom.name.strip() == args.patom :
            memz.append(atom.z)
    memcent = np.asarray(memz).mean()
    print "Membrane center = %.3f"%memcent

    water_z = []
    water_resid = []
    for i, residue in enumerate(struct.residues) :
        if residue.atoms[0].name.strip() == args.watom :
            water_z.append(residue.atoms[0].z)
            water_resid.append(i)

    water_z = np.asarray(water_z)
    sortidx = np.argsort(water_z - memcent)
    remove_these = [water_resid[i] for i in sortidx[:args.n[0]]]
    remove_these2 = [water_resid[i] for i in sortidx[-args.n[1]:]]
    remove_these.extend(remove_these2)

    residues = []
    atoms = []
    for i, residue in enumerate(struct.residues) :
        if i not in remove_these :
            residues.append(residue)
            atoms.extend(residue.atoms)
    struct.atoms = atoms
    struct.residues = residues

    if not args.nobox :
        zcoord = np.asarray([atom.z for atom in atoms])
        zmin = np.floor(zcoord.min()) - 1.0
        zlen = np.floor(zcoord.max() - zcoord.min() + 2.0)
        for atom in atoms :
            atom.z -= zmin
            atom.xyz[2] -= zmin
        struct.box[2] = zlen

    struct.write(args.out)
