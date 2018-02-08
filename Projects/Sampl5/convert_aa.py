# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Scripts used to convert the provided all-atom system input in Lammps format
to a hybrid description. The solute is retained as all-atom, whereas the
solvent molecules are made into CG.

The script creates the following files:
1) A PDB-file for the solute
2) A Lammps datafile for the solute, include force field parameters
3) A Lammps datafile for the CG solvent box

The script assumes that it is given a Lammps datafile and parses the
solute name and solvent from the filename. It also assumes that forcefield
parameters can be found in a similarly named file but with the "input" extension.

"""

import sys
import os

import numpy as np

from sgenlib import lammps
from sgenlib import pdb
from sgenlib import geo

def _make_solutepdb(atoms, masses, ligand):
    """
    Make a solute PDB file
    """
    mass2elem = {1.0080000:"H",
                 12.0100000:"C",
                 14.0100000:"N",
                 16.0000000:"O",
                 32.0600000:"S",
                 35.4500000:"Cl",
                 19.0000000:"F",
                 79.9000000:"Br"}

    pdbout = pdb.PDBFile()
    res = pdb.Residue()
    res.resname = ligand
    res.serial = 1
    for i, (a, m) in enumerate(zip(atoms, masses),1):
        atom = pdb.Atom()
        atom.serial = i
        atom.hetatm = False
        atom.name = "%s%d"%(mass2elem[m],i) # Atom names from masses and serial number
        atom.resname = ligand
        atom.residue = 1
        atom.set_xyz(a.xyz)
        res.atoms.append(atom)
        pdbout.atoms.append(atom)
    pdbout.residues.append(res)
    pdbout.write(ligand+".pdb")

def _make_solutedata(aadata, solutedata, solatype, solidx, pairparams,
                        ligand, infilename):
    """
    Make a solute datafile, by copy everything from the all-atom
    datafile that belongs to the solute
    """
    # Check the bonds, copying over bonds in solute and record the bond type
    solbondtypes = []
    for bond in aadata.bonds :
        onsolute = False
        for a in bond.atoms :
            if a in solidx:
                onsolute = True
                break
        if onsolute :
            solutedata.bonds.append(bond)
            solbondtypes.append(bond.param)

    # Same for angles...
    solangletypes = []
    for angle in aadata.angles :
        onsolute = False
        for a in angle.atoms :
            if a in solidx:
                onsolute = True
                break
        if onsolute :
            solutedata.angles.append(angle)
            solangletypes.append(angle.param)

    # Same for dihedrals...
    soldihedtypes = []
    for dihedral in aadata.dihedrals :
        onsolute = False
        for a in dihedral.atoms :
            if a in solidx:
                onsolute = True
                break
        if onsolute :
            solutedata.dihedrals.append(dihedral)
            soldihedtypes.append(dihedral.param)

    # Make new atom types
    solutedata.atomtypes = [t for t in aadata.atomtypes if t.idx in solatypes]
    for t in solutedata.atomtypes:
        t.epsilon = pairparams[t.idx].epsilon
        t.sigma = pairparams[t.idx].sigma

    # Make new bond, angle and dihedral types
    solutedata.bondtypes = [t for t in aadata.bondtypes if t.idx in solbondtypes]
    for t in solutedata.bondtypes:
        t.func = "" # This if to work with insert_in_elba.py, which assumes harmonic function
    solutedata.angletypes = [t for t in aadata.angletypes if t.idx in solangletypes]
    solutedata.dihedraltypes = [t for t in aadata.dihedraltypes if t.idx in soldihedtypes]

    solutedata.title = "Made with convert_aa.py from %s"%os.path.basename(infilename)
    solutedata.box = np.array(aadata.box,copy=True)
    solutedata.write("data.%s"%ligand,writeparams=True)

def _make_waterdata(atoms, box, ligands, infilename):
    """
    Make a datafile for a CG box with water molecules
    """
    outdata = lammps.Datafile()
    atoms = atoms[::3] # Map to the oxygen

    for i, aaatom in enumerate(atoms, 1):
        cgatom = lammps.Atom()
        cgatom.idx = i
        cgatom.atype = 1
        cgatom.set_xyz(aaatom.xyz) # This is the only thing used from all-atom
        cgatom.charge = 0
        cgatom.set_mu(geo.sphere_rand(0.541)) # Random dipole direction
        cgatom.diameter = 3.0
        cgatom.density = 2.7
        cgatom.kind = "cg"
        outdata.atoms.append(cgatom)

    outdata.atomtypes = [None]*6
    outdata.bondtypes = [None]*5
    outdata.angletypes = [None]*5

    outdata.box = np.array(box, copy=True)
    outdata.title = "Converted to CG with convert_aa.py from %s"%os.path.basename(infilename)
    outdata.write("data.%s-waterbox"%ligand)

def _make_cyclohexanedata(atoms, box, ligands, infilename):
    """
    Make a datafile for a CG box with cyclohexane molecules
    """

    outdata = lammps.Datafile()
    c1_atoms = atoms[::18] # Map to first carbon
    c2_atoms = atoms[1::18] # Map to second carbon
    c3_atoms = atoms[3::18] # Map to fourth carbon

    for i, (c1atom, c2atom, c3atom) in enumerate(zip(c1_atoms,c2_atoms,c3_atoms)):
        aatoms = [c1atom,c2atom,c3atom]
        for j in range(3):
            cgatom = lammps.Atom()
            cgatom.idx = i*3+j+1
            cgatom.atype = 6
            cgatom.set_xyz(aatoms[j].xyz)
            cgatom.charge = 0
            cgatom.set_mu([0,0,0])
            cgatom.diameter = 4.5
            cgatom.density = 0.9
            cgatom.kind = "cg"
            cgatom.molecule = i+1
            outdata.atoms.append(cgatom)
        outdata.bonds.append(lammps.Connectivity("%d 5 %d %d"%(len(outdata.bonds)+1,i*3+1,i*3+2)))
        outdata.bonds.append(lammps.Connectivity("%d 5 %d %d"%(len(outdata.bonds)+1,i*3+1,i*3+3)))
        outdata.bonds.append(lammps.Connectivity("%d 5 %d %d"%(len(outdata.bonds)+1,i*3+3,i*3+2)))

    outdata.atomtypes = [None]*6
    outdata.bondtypes = [None]*5
    outdata.angletypes = [None]*5

    outdata.box = np.array(box, copy=True)
    outdata.title = "Converted to CG with convert_aa.py from %s"%os.path.basename(infilename)
    outdata.write("data.%s-chexanebox"%ligand)

if __name__ == '__main__' :

    filename = sys.argv[1]
    filehead = os.path.splitext(os.path.basename(filename))[0]
    solvent = filehead.split("-")[0]
    ligand = filehead.split("_")[1]

    aadata = lammps.Datafile(filename)
    aapairparams = lammps.Includefile(os.path.splitext(filename)[0]+".input")

    solutedata = lammps.Datafile()
    solatypes = []
    solidx = []
    solutemasses = []
    solventatoms = []
    # Split up solute and solvent atoms
    for atom in aadata.atoms:
        if atom.molecule == 1 :
            solatypes.append(atom.atype)
            solidx.append(atom.idx)
            solutedata.atoms.append(atom)
            solutemasses.append(aadata.atomtypes[atom.atype-1].mass)
        else:
            solventatoms.append(atom)
    solatypes = set(solatypes)

    _make_solutepdb(solutedata.atoms, solutemasses, ligand)
    _make_solutedata(aadata, solutedata, solatypes, solidx,
                        aapairparams.lj_same_dict(), ligand, filename)

    if solvent == "water":
        _make_waterdata(solventatoms, aadata.box, ligand, filename)
    elif solvent == "cyclohexane":
        _make_cyclohexanedata(solventatoms, aadata.box, ligand, filename)
