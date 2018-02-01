# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to build a CG or AA/CG dual-resolution
system from an AA system. Converts AA residues
automatically to CG using dictionary of conversions.

Writes out LAMMPS datafiles and include file
for the force field, as well as a PDB-file.
"""

import sys
import math
import random
import argparse
import os
import copy

import numpy as np
from scipy.spatial.distance import cdist

from sgenlib import lammps
from sgenlib import pdb


def _ntypes(array):
    """
    Find the number of types from a list (could be masses or connectivities)
    """
    n = 0
    for m in array:
        n = max(n, m.idx)
    return n


def _generate_aa_residue(residue, molidx, resdata, sysdata):
    """
    Generates an aa residue by copying most of the structure
    from a datafile template, but the coordinates from a PDB residue
    """
    n = len(sysdata.atoms)
    for i, (ratom, datom) in enumerate(zip(residue.atoms, resdata.atoms)):
        atom = copy.deepcopy(datom)
        atom.idx = atom.idx + n
        atom.set_xyz(ratom.xyz)
        atom.molecule = molidx
        sysdata.atoms.append(atom)

    for bond in resdata.bonds:
        b = lammps.Connectivity(record="%d %d %d %d" % (len(sysdata.bonds) + 1, bond[0], bond[1] + n, bond[2] + n))
        sysdata.bonds.append(b)

    for angle in resdata.angles:
        a = lammps.Connectivity(record="%d %d %d %d %d" % (len(sysdata.angles) + 1, angle[0], angle[1] + n, angle[2] + n, angle[3] + n))
        sysdata.angles.append(a)

    for dihedral in resdata.dihedrals:
        a = lammps.Connectivity(record="%d %d %d %d %d %d" % (len(sysdata.dihedrals) + 1, dihedral[0], dihedral[1] + n, dihedral[2] + n, dihedral[3] + n, dihedral[4] + n))
        sysdata.dihedrals.append(a)

def _generate_aa_wat(residue, molidx, sysdata, include):

    n = len(sysdata.atoms)
    charges = [ -0.834, 0.417, 0.417]
    mass = [9.9514, 3.024, 3.024]
    nt = len(include.masses)
    types = [nt-1,nt,nt]
    for i, ratom in enumerate(residue.atoms):
        atom = lammps.Atom()
        atom.idx = n + i + 1
        atom.set_xyz(ratom.xyz)
        atom.molecule = molidx
        atom.q = charges[i]
        atom.atype = types[i]
        atom.comment = "# water"
        atom.diameter = 0.0
        atom.density = mass[i]
        atom.set_mu([0.0, 0.0, 0.0])
        sysdata.atoms.append(atom)

    nb = len(include.bondparams)
    b = lammps.Connectivity(record="%d %d %d %d" % (len(sysdata.bonds) + 1, nb, n+1, n+2))
    b.comment = "# O-H bond"
    sysdata.bonds.append(b)
    b = lammps.Connectivity(record="%d %d %d %d" % (len(sysdata.bonds) + 1, nb, n+1, n+3))
    b.comment = "# O-H bond"
    sysdata.bonds.append(b)

    na = len(include.angleparams)
    a = lammps.Connectivity(record="%d %d %d %d %d" % (len(sysdata.angles) + 1, na, n+2, n+1, n+3))
    a.comment = "# H-O-H angle"
    sysdata.angles.append(a)


def _add_water_param(datafile):

    n = _ntypes(datafile.atomtypes)
    # Oxygen type
    datafile.atomtypes.append(lammps.AtomType())
    datafile.atomtypes[-1].idx = n + 1
    datafile.atomtypes[-1].mass = 9.9514
    datafile.atomtypes[-1].epsilon = 0.1521
    datafile.atomtypes[-1].sigma = 3.1507
    # Hydrogen type
    datafile.atomtypes.append(lammps.AtomType())
    datafile.atomtypes[-1].idx = n + 2
    datafile.atomtypes[-1].mass = 3.024
    datafile.atomtypes[-1].epsilon = 0.0
    datafile.atomtypes[-1].sigma = 0.0

    # O-H bond type
    n = _ntypes(datafile.bondtypes)
    datafile.bondtypes.append(lammps.ConnectivityParam())
    datafile.bondtypes[-1].idx = n + 1
    datafile.bondtypes[-1].params = ["450.0", "0.9572"]

    # H-O-H angle type
    n = _ntypes(datafile.angletypes)
    datafile.angletypes.append(lammps.ConnectivityParam())
    datafile.angletypes[-1].idx = n + 1
    datafile.angletypes[-1].params = ["55.0", "104.52"]

if __name__ == '__main__':

    # Command-line input
    parser = argparse.ArgumentParser(description="Converting a AA to CG or AA/CG")
    parser.add_argument('file', help="the PDB or GRO file")
    parser.add_argument('-i', '--include', help="the LAMMPS include file")
    parser.add_argument('-o', '--out', help="the output prefix", default="converted")
    parser.add_argument('-b', '--box', type=float, nargs="+", help="the box dimensions", default=[0.0, 0.0, 0.0])
    parser.add_argument('-a', '--atomistic', nargs="+", help="data file(s) for atomistic solutes", default=[])
    parser.add_argument('-c', '--converter', help="the dictionary with conversion rules")
    parser.add_argument('-p', '--pairfunc', help="the pair function for the AA", default="lj/charmm/coul/long")
    parser.add_argument('-w', '--watrad', type=float, help="the water radius to keep atomistic")
    args = parser.parse_args()

    # Load a converter
    converter = lammps.Aa2Cg()
    if args.converter is None:
        converter.read(lammps.get_filename("aa2cg.dat"))  # The default
    else:
        converter.read(args.converter)

    # Create a Datafile and PDBFile
    pdbfile = pdb.PDBFile(args.file)  # Input PDB
    takeres = [True for res in pdbfile.residues]
    data = lammps.Datafile()
    pdbout = pdb.PDBFile()  # Output PDB

    # Load the force field file
    include = lammps.Includefile(args.include)

    if args.atomistic:
        # At the moment, multiple AA solutes are not supported with atomistic water radius
        if len(args.atomistic) > 1 :
            args.watrad = None
        natomtypes = _ntypes(include.masses)
        contypes = [_ntypes(include.bondparams), _ntypes(include.angleparams), _ntypes(include.dihedralparams)]
        # Set the functional form of the CG particles if will retain some atomistic molecules
        for pair in include.pair_coeff:
            pair.func = "lj/sf/dipole/sf"
            pair.hybrid = 2

    # Load datafiles for given solutes, will assume these are atomistic
    aa_datafiles = {}
    aa_range = None
    for sol in args.atomistic:
        res, filename = sol.split("=")
        res = res.lower()
        if res[0] == ":":
            aa_range = res
        aa_datafiles[res] = lammps.Datafile(filename)
        if args.watrad is not None:
            _add_water_param(aa_datafiles[res])
        # Extend the force field parameters
        # by extending the inclusion file, we will automatically update the parameter index
        ELBA_FUNC = "lj/sf/dipole/sf"
        include.extend_from_data(aa_datafiles[res], lj_hybrid=-1, lj_func=args.pairfunc,
                                 lj_hybrid_mix={ELBA_FUNC: 1, args.pairfunc: 2}, lj_func_mix={ELBA_FUNC: ELBA_FUNC, args.pairfunc: args.pairfunc}, ang_func="harmonic")
        # Update the atom and conectivity parameters
        for atom in aa_datafiles[res].atoms:
            atom.atype = atom.atype + natomtypes
            atom.diameter = 0.0
            atom.density = aa_datafiles[res].atomtypes[atom.atype - natomtypes - 1].mass
            atom.set_mu([0.0, 0.0, 0.0])
        conlist = [aa_datafiles[res].bonds, aa_datafiles[res].angles, aa_datafiles[res].dihedrals]
        for cons, ntypes in zip(conlist, contypes):
            for con in cons:
                con.param = con.param + ntypes
        natomtypes = _ntypes(include.masses)

    # Add a range of all-atom residues to the data file
    moli = 0
    if aa_range is not None:
        first, last = map(lambda x: int(x) - 1, aa_range[1:].split("-"))
        resall = pdb.Residue()
        allres = []
        for i in range(first, last + 1):
            for atom in pdbfile.residues[i].atoms:
                resall.append(atom)
            allres.append(pdbfile.residues[i])
            takeres[i] = False
        _generate_aa_residue(resall, 1, aa_datafiles[aa_range], data)
        moli = 1

        if args.watrad is not None:
            ninside = 0
            com = [np.asarray([res.collect("centerofmass") for res in allres]).mean(axis=0)]
            watres = []
            mindist = []
            for i,residue in enumerate(pdbfile.residues) :
                if not residue.resname in ["HOH","WAT","SOL"] : continue
                mindist.append(np.sum((com-residue.atoms[0].xyz)**2))
                #xyz1 = pdbfile.xyz[residue.atoms[0].idx:residue.atoms[-1].idx+1,:]
                #print xyz1[0],
                #mindist.append(cdist(xyz1,com,"sqeuclidean").min())
                watres.append((i,residue))
            mindist = np.asarray(mindist)
            for ri in np.argsort(mindist)[:args.watrad]:
                i,residue = watres[ri]
                _generate_aa_wat(residue, 2, data, include)
                allres.append(residue)
                takeres[i] = False
                ninside += 1
            moli = 2
            print "Kept %d water molecules as all-atom"%ninside

        pdbout.extend_residues(allres, dochains=False)

    # Convert residues
    all_coords = []
    nwat = 0
    for i, (res, takethis) in enumerate(zip(pdbfile.residues, takeres)):
        if not takethis:
            continue
        moli += 1
        res2 = res.resname.strip().lower()
        found = False
        # If we have an all-atom datafile as a template, keep it as all-atom
        if res2 in aa_datafiles:
            _generate_aa_residue(res, moli - nwat, aa_datafiles[res2], data)
            coord = res.collect("xyz")
            pdb.make_pdbres(coord, [atom.name for atom in res.atoms], res2, pdbout)
            found = True
        # Otherwise convert it to CG
        else:
            for residue in converter.residues:
                if residue.name == res2:
                    coord = residue.generate_cg(res, moli, data)
                    all_coords.extend(coord)
                    pdb.make_pdbres(coord, residue.cg_names, res2, pdbout)
                    found = True
                    break
        # If we could not find a conversion, we will convert the residue to a water bead
        if not found:
            for residue in converter.residues:
                if residue.name == "wat":
                    nwat = nwat + 1
                    coord = residue.generate_cg(res, 0, data)
                    all_coords.extend(coord)
                    pdb.make_pdbres(coord, residue.cg_names, "wat", pdbout)

    all_coords = np.array(all_coords)

    print "Minimum of coordinates = %.3f %.3f %.3f" % tuple(all_coords.min(axis=0))
    print "Maximum of coordinates = %.3f %.3f %.3f" % tuple(all_coords.max(axis=0))
    print "Average of coordinates = %.3f %.3f %.3f" % tuple(all_coords.mean(axis=0))

    # Settings the correct number of atom and connectivity types
    data.atomtypes = [None] * len(include.masses)
    data.bondtypes = [None] * len(include.bondparams)
    data.angletypes = [None] * len(include.angleparams)
    data.dihedraltypes = [None] * len(include.dihedralparams)

    # Setting the box of the datafile
    if len(args.box) == 6 :
        data.box = args.box
    else :
        if all_coords.mean(axis=0).sum() > 10:  # Checking if center is at origin or not
            data.box = [0.0, 0.0, 0.0, args.box[0], args.box[1], args.box[2]]
        else:
            data.box = [-args.box[0] / 2.0, -args.box[1] / 2.0, -args.box[2] / 2.0, args.box[0] / 2.0, args.box[1] / 2.0, args.box[2] / 2.0]

    # Setting correct type for all atoms
    if args.atomistic:
        for atom in data.atoms:
            atom.kind = "cg/aa"

    # Setting box and adding connectivity to PDB-file
    pdbout.box = args.box

    def add_con(pdbfile, f):
        for bnd in data.bonds:
            f.write("CONECT%5d%5d\n" % (bnd.atoms[0], bnd.atoms[1]))

    # Write out datafile, pdbfile and force field
    print "Saving LAMMPS data file to %s" % ("data." + args.out)
    data.write("data." + args.out)
    print "Saving PDB file to %s" % (args.out + ".pdb")
    pdbout.write(args.out + ".pdb", add_extra=add_con)
    if args.atomistic:
        print "Saving LAMMPS inclusion file to %s" % ("forcefield." + args.out)
        include.write("forcefield." + args.out)
