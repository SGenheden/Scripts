# Author: Samuel Genheden samuel.genheden@gmail.com
"""
Program to extend the structure of a lipid such that
the acyl chains are longer. Assumes Slipid/Charmm names

Adds the new atoms to idealized coordinates, assuming sp3 hybridization

Uses in membrane engineering project

Examples
--------
extend_chain.py -f mem.gro -a 20 -r DOPC DOPC DOPC -n 16 8 4 -w AOPC BOPC LOPC -a 20 22 24 -c sn1
    extends the sn1 chain of DOPC residues, 16 will be 20 in length and called AOPC,
    8 will be 22 in length and called BOPC and 4 will be 24 in length and called LOPC
"""

import sys
import argparse

import numpy as np
import numpy.random as random

from sgenlib import pdb
from sgenlib import geo

def _create_atoms(residue, allatoms, oldend, template, htemplate) :
    """
    Add new atoms to a residue and change a hydrogen to a carbon

    Parameters
    ----------
    residue : pdb.Residue
        the residue to extend
    allatoms : list of pdb.Atom
        the list of all atoms
    oldend : integer
        the current terminal carbon
    template : string
        the name prefix of the carbon atoms
    htemplate: string
        the name postfox of the hydrogen atoms
    """
    # Find the preceeding carbons
    atoms = []
    for idx in range(oldend-2, oldend+1) :
        for atom in residue.atoms :
            if atom.name == "%s%d"%(template,idx) :
                atoms.append(atom)

    # Find the hydrogen atom and change it to a carbon
    hatoms = [None]*3
    angles = [0]*3
    hatom = None
    for atom in residue.atoms :
        for i, hi in enumerate(htemplate) :
            if atom.name == "H%d%s"%(oldend, hi) :
                hatoms[i] = atom
                angles = np.asarray([0,0,1.0])*(atoms[2].xyz-atom.xyz)
    hatom = hatoms[2]
    hatom.name = "%s%d"%(template, oldend+1)
    # Update the coordinates of the new terminal carbon
    dihed = geo.dihedral_protoms(hatom.xyz, atoms[2].xyz, atoms[1].xyz, atoms[0].xyz)*180.0/np.pi
    hatom.set_xyz(geo.build_xyz(atoms[2].xyz, atoms[1].xyz, atoms[0].xyz,
                    1.53, 114.0, dihed))

    # Now create 3 new hydrogen atoms
    for i in range(3) :
        new_atom = pdb.Atom(record=hatom.__str__())
        new_atom.name = "H%d%s"%(oldend+1,htemplate[i])
        # Use idealized coordinates
        new_atom.set_xyz(geo.build_xyz(hatom.xyz, atoms[2].xyz, atoms[1].xyz,
                        1.1, 110.0, -60+120*i))
        allatoms.insert(allatoms.index(hatom)+1+i, new_atom)
        residue.atoms.append(new_atom)

def _extend_res(residue, atoms, lastatom, chains, resname) :
    """
    Extend the structure of a lipid molecule

    Parameters
    ----------
    residue : pdb.Residue
        the residue to extend
    atoms : list of pdb.Atom
        the list of all atoms
    lastatom : integer
        extend the acyl chains to this length
    chains : list of strings
        the chains to extend, sn1 or sn2 or both
    resname : string
        the new residue name
    """

    # Find the current terminal carbon for the sn-1 and sn-2 chains
    found2 = 1
    found3 = 1
    for idx in range(2,20) :
        for atom in residue.atoms :
            if atom.name == "C2%d"%idx :
                found2 = idx
            elif atom.name == "C3%d"%idx :
                found3 = idx

    # ... and then extend them
    if "sn2" in chains :
        for add_idx in range(found2, lastatom) :
            _create_atoms(residue, atoms, add_idx, "C2", ["R", "S", "T"])
    if "sn1" in chains :
        for add_idx in range(found3, lastatom) :
            _create_atoms(residue, atoms, add_idx, "C3", ["X", "Y", "Z"])

    for atom in residue.atoms :
        atom.resname = resname

if __name__ == '__main__':

    # We will displace the upper part of the membrane this much to accomodate the longer chains
    Z_CONST = np.asarray([0.0, 0.0, 6])

    parser = argparse.ArgumentParser(description="Program to extend the acyl chains of lipids")
    parser.add_argument('-f','--file',help="the input gro-file")
    parser.add_argument('-r','--resname', nargs="+",help="the residue names of the lipids to extend")
    parser.add_argument('-n','--nreplace',nargs="+",type=int,help="how many lipids molecules to replace")
    parser.add_argument('-w','--newname',nargs="+",help="the new residue name")
    parser.add_argument('-a','--atom',nargs="+",type=int,help="the new end atom")
    parser.add_argument('-c','--chains',nargs="+",help="the chains to extend",default=["sn1", "sn2"])
    parser.add_argument('--seed', type=int, help="the random seed", default=1124142)
    parser.add_argument('-o', '--out', help="the output name", default="extended.gro")
    args = parser.parse_args()

    random.seed(args.seed)

    struct = pdb.PDBFile(args.file)

    # Displace all atoms in the upper leaflet
    midz = struct.xyz[:,2].mean()
    for residue in struct.residues :
        if residue.collect("centerofmass")[2] > midz :
            for atom in residue.atoms :
                atom.set_xyz(atom.xyz+Z_CONST)
    struct.box += Z_CONST

    # Find residues in lower and upper part of the membrane
    res_low = []
    res_upp = []
    for residue in struct.residues :
        if residue.collect("centerofmass")[2] > midz :
            res_upp.append(residue)
        else :
            res_low.append(residue)
        residue.taken = False


    for iresname, inreplace, iatom, inewname in zip(args.resname, args.nreplace, args.atom, args.newname) :

        ires_low = [res for res in res_low if res.resname == iresname and not res.taken]
        ires_upp = [res for res in res_upp if res.resname == iresname and not res.taken]
        ntaken = 0

        # Extend nreplace/2 number of molecules in the lower part of the membrane
        while ntaken < int(np.ceil(inreplace * 0.5)) :
            # Find a residue that is not taken
            takei = random.randint(0,len(ires_low))
            while ires_low[takei].taken : takei = random.randint(0,len(ires_low))
            _extend_res(ires_low[takei], struct.atoms, iatom, args.chains, inewname)
            ires_low[takei].taken = True
            ntaken += 1

        # and now the upper part of the membrane
        while ntaken < inreplace :
            takei = random.randint(0,len(ires_upp))
            while ires_upp[takei].taken : takei = random.randint(0,len(ires_upp))
            _extend_res(ires_upp[takei], struct.atoms, iatom, args.chains, inewname)
            ires_upp[takei].taken = True
            ntaken += 1

        print "Replaced %d %s residues"%(inreplace, iresname)
    struct.write(args.out)
