# Author: Samuel Genheden samuel.genheden@gmail.com
"""
Program to replace double bond in acyl chain with single bond

Uses in membrane engineering project

saturate_chain.py -f extended.gro -r AOPC BOPC LOPC
    replaces the unsaturated sn-1 chain in AOPC BOPC and LOPC with saturated chains
"""

import argparse

import numpy as np

from sgenlib import geo
from sgenlib import pdb

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Program to saturate chains")
    parser.add_argument('-f','--file',help="the input gro-file")
    parser.add_argument('-r','--residues',nargs="+",help="the residues to replace")
    parser.add_argument('-o', '--out', help="the output name", default="saturated.gro")
    args = parser.parse_args()

    struct = pdb.PDBFile(args.file)

    selres = [res for res in struct.residues if res.resname in args.residues]

    for i, residue in enumerate(selres) :
        carbons = []
        hydrogens = []
        for atom in residue.atoms :
            if atom.name in ["C37", "C38", "C39", "C310","C311"] :
                carbons.append(atom)
            elif atom.name in ["H9X", "H10X"] :
                hydrogens.append(atom)

        # This is the dihedral to the next carbon atom after the hydrogen
        t1 = geo.dihedral_protoms(carbons[0].xyz, carbons[1].xyz, carbons[2].xyz, carbons[3].xyz)*180.0/np.pi
        t2 = geo.dihedral_protoms(carbons[1].xyz, carbons[2].xyz, carbons[3].xyz, carbons[4].xyz)*180.0/np.pi

        # Move the hydrogen in the unsaturated chain
        hydrogens[0].set_xyz(geo.build_xyz(carbons[2].xyz, carbons[1].xyz, carbons[0].xyz,
                        1.1, 110.0, t1+120))
        hydrogens[1].set_xyz(geo.build_xyz(carbons[3].xyz, carbons[2].xyz, carbons[1].xyz,
                        1.1, 110.0, t2+120))

        # Create two new atoms, making the chain saturated
        new = pdb.Atom(record=hydrogens[0].__str__())
        new.set_xyz(geo.build_xyz(carbons[2].xyz, carbons[1].xyz, carbons[0].xyz,
                        1.1, 110.0, t1-120))
        new.name = "H9Y"
        struct.atoms.insert(struct.atoms.index(hydrogens[0])+1, new)
        residue.atoms.append(new)
        new.resname = residue.resname

        new = pdb.Atom(record=hydrogens[1].__str__())
        new.set_xyz(geo.build_xyz(carbons[3].xyz, carbons[2].xyz, carbons[1].xyz,
                        1.1, 110.0, t2-120))
        new.name = "H10Y"
        struct.atoms.insert(struct.atoms.index(hydrogens[1])+1, new)
        residue.atoms.append(new)
        new.resname = residue.resname

    struct.renumber(doatoms=True, doresidues=False)
    struct.write(args.out)
