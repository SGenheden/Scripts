# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to insert cholesterol in a CG membrane by replacing phospholipids

The current mapping is:
    C2 --> C3B
    C1 --> C2B
    R5 --> C1B
    R2 --> C1A
    R1 --> GL2
    ROH --> GL1

Example:
    python insert_chol.py -f popc_mem.gro -c chol_single.gro -n 60
"""

import argparse
import copy

import numpy as np

from sgenlib import pdb
from sgenlib import fitting

def _do_fit(cholstruct, lipidstruct) :

    fit_map = [["C2","C3B"],["C1","C2B"],["R5","C1B"],["R2","C1A"],["R1","GL2"],["ROH","GL1"]]
    ref = np.zeros([len(fit_map),3])
    mob = np.zeros([len(fit_map),3])
    for i, pair in enumerate(fit_map) :
        for atom in cholstruct.atoms :
            if atom.name.strip() == pair[0]:
                mob[i,:] = atom.xyz
                break
        for atom in lipidstruct.atoms :
            if atom.name.strip() == pair[1]:
                ref[i,:] = atom.xyz
                break
    xyz = cholstruct.collect("xyz")

    xyz = fitting.dofit(ref, mob, xyz)
    cholstruct2 = copy.deepcopy(cholstruct)
    cholstruct2.update_xyz(xyz)

    lipidstruct.resname = cholstruct2.resname
    lipidstruct.atoms = cholstruct2.atoms

def _extract_lipids(residues, lipidnames) :

    phosphates = []
    lipids = []
    for residue in residues :
        if residue.resname in lipidnames :
            patom = None
            for atom in residue.atoms:
                if atom.name.strip() == "PO4":
                    patom = atom
                    break
            if patom is not None:
                phosphates.append(patom)
                lipids.append(residue)
    return lipids, phosphates

def _replace(lipids, phosphates, cholstruct, nreplace, taken, leaffunc) :

    ntaken = 0
    while ntaken < nreplace :
        probe = np.random.randint(0, len(lipids))
        while taken[probe] or leaffunc(phosphates[probe]) :
            probe = np.random.randint(0, len(lipids))
        ntaken = ntaken + 1
        taken[probe] = True
        _do_fit(cholstruct, lipids[probe])

if __name__ == '__main__' :

    parser = argparse.ArgumentParser(description="Insert cholesterol by replacing lipids")
    parser.add_argument('-f','--file',help="the lipid structure file")
    parser.add_argument('-c','--chol',help="the cholesterol template structure")
    parser.add_argument('-o','--out',help="the output file",default="chol.gro")
    parser.add_argument('-l','--lipids', nargs="+", help="the lipid residue names", default=["POPC"])
    parser.add_argument('-n','--nreplace', type=int, help="the number of lipids to replace")
    args = parser.parse_args()

    cholstruct = pdb.PDBFile(args.chol).residues[0]
    memstruct = pdb.PDBFile(args.file)
    lipids, phosphates = _extract_lipids(memstruct.residues, args.lipids)

    midz = np.asarray([p.z for p in phosphates]).mean()
    taken = [False]*len(lipids)
    nupper = int(0.5*args.nreplace)

    _replace(lipids, phosphates, cholstruct, nupper, taken, lambda p: p.z < midz)
    _replace(lipids, phosphates, cholstruct, args.nreplace - nupper, taken, lambda p: p.z > midz)
    cholesterols = [lipid for ltaken, lipid in zip(taken, lipids) if ltaken]

    # Now we need to fix the structure a litle bit, sort it etc.

    # First find out where to insert the cholesterol, right before the first water
    wati = 0
    for i, res in enumerate(memstruct.residues) :
        if res.resname.strip() == "W":
            wati = i
            break

    # Remove the cholesterols from their old lipid position to just before the waters
    for cholesterol in cholesterols :
        memstruct.residues.remove(cholesterol)
        memstruct.residues.insert(wati-1, cholesterol)

    # Re-insert cholesterol atoms and renumber atoms and residues
    memstruct.atoms = []
    for ri, residue in enumerate(memstruct.residues, 1):
        for atom in residue.atoms:
            atom.serial = len(memstruct.atoms)+1
            atom.residue = ri
            memstruct.atoms.append(atom)

    memstruct.write(args.out)
