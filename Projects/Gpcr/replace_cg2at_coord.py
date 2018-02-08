# Author: Samuel Genheden samuel.genheden@gmail.com

import sys
import os

from sgenlib import pdb


def atom_by_name(residue, atom) :
    atom = atom.strip()
    if atom[0] in ["1","2","3"]:
        atom = atom[1:]+atom[0]
    if residue.resname == "ILE" and atom == "CD" :
        atom = "CD1"
    return residue.atom_by_name(atom)

converted_pdb = pdb.PDBFile(sys.argv[1])
crystal_pdb = pdb.PDBFile(sys.argv[2])

nprotres = len(crystal_pdb.residues)

for conv_res, crys_res in zip(converted_pdb.residues[:nprotres], crystal_pdb.residues):
    for conv_atom in conv_res.atoms:
        crys_atom = atom_by_name(crys_res, conv_atom.name)
        if crys_atom is None:
            print "Could not find %s in %s%d"%(conv_atom.name, conv_atom.resname, conv_atom.residue)
        else:
            conv_atom.set_xyz(crys_atom.xyz)

name, ext = os.path.splitext(sys.argv[1])
converted_pdb.write(name+"-crystal"+ext)
