# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to make PQR file from a Gromacs topology and a structure file
"""

import argparse
import os
import sys

from sgenlib import pdb
from sgenlib import gmx

def _parse_radii(atom, topol) :

    def _neigh(iatom, topol, hydrogen=False):
        neigh = []
        for bond in topol.bonds :
            if len(bond.atoms) != 2 : continue
            if bond.atoms[0] == iatom  :
                neigh.append(bond.atoms[1])
                if hydrogen : break
            elif bond.atoms[1] == iatom :
                neigh.append(bond.atoms[0])
                if hydrogen : break
        return neigh

    if atom.mass == 1.008 :
        hneigh = _neigh(atom.id, topol, hydrogen=True)
        atom_neigh = topol.atoms[hneigh[0]-1]
        if atom_neigh.mass in [12.01, 14.01]:
            if len(_neigh(atom_neigh.id, topol)) == 4 :
                return 0.0
            else :
                return 1.0
        else :
            return 1.0
    elif atom.mass == 12.01 :
        if len(_neigh(atom.id, topol)) == 4 :
            return 2.0
        else :
            return 1.7
    elif atom.mass == 14.01 :
        if len(_neigh(atom.id, topol)) == 4 :
            return 2.0
        else :
            return 1.5
    elif atom.mass == 16.00 :
        return 1.4
    elif atom.mass == 32.06 :
        return 1.85
    else :
        return 2.0

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Program make a PQR file")
    parser.add_argument('-f','--file',help="the input gro-file")
    parser.add_argument('-p','--topol',help="the topology file")
    parser.add_argument('-o','--out',help="the output pqr-file",default="conv.pqr")
    parser.add_argument('-m','--mol',nargs="+",help="the molecule to write out",default=["protein"])
    args = parser.parse_args()

    ref = gmx.TopFile(args.topol)
    pdbfile = pdb.PDBFile(args.file)

    mols = []
    for topol_mol in ref.moleculetypes :
        for argmol in args.mol :
            if topol_mol.name.lower() == argmol.lower() :
                mols.append(topol_mol)
                break

    natom = sum([len(mol.atoms) for mol in mols])
    pdbfile.residues = [res for res in pdbfile.residues if res.atoms[0].serial <= natom]
    pdbfile.atoms = pdbfile.atoms[:natom]

    start = 0
    for mol in mols  :
        patoms = pdbfile.atoms[start:start+len(mol.atoms)]
        for patom, tatom in zip(patoms, mol.atoms):
            patom.occupancy = tatom.charge
            patom.bfactor = _parse_radii(tatom, mol)
            patom.pqratom = True
        start = len(mol.atoms)
    pdbfile.write(args.out)
