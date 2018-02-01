# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to truncate a data file to a QM system
Part of ComQum-ELBA
"""

import argparse
from collections import namedtuple

import numpy as np

from sgenlib import lammps

def _trunc_connectivity(conlist, atom_ids) :

    newlist = []
    for con in conlist:
        keep = False
        for atom in con.atoms :
            if atom in atom_ids :
                keep = True
                break
        if keep :
            con.atoms = [atom_ids.index(atom)+1 for atom in con.atoms]
            con.idx = len(newlist)+1
            newlist.append(con)
    return newlist

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Truncate a LAMMPS datafile")
    parser.add_argument('file', help="the data file.")
    parser.add_argument('-m', '--mol', type=int, help="the molecule id of the QM system")
    parser.add_argument('-o', '--out',help="the output file",default="data.sys1")
    args = parser.parse_args()

    if args.file is None:
        print "No input file specified. Exiting!"
        quit()

    datafile = lammps.Datafile(args.file)
    new_atoms = [atom for atom in datafile.atoms if atom.molecule == args.mol]
    atom_ids = [atom.idx for atom in new_atoms]

    datafile.bonds = _trunc_connectivity(datafile.bonds, atom_ids)
    datafile.angles = _trunc_connectivity(datafile.angles, atom_ids)
    datafile.dihedrals = _trunc_connectivity(datafile.dihedrals, atom_ids)

    for i, atom in enumerate(new_atoms, 1):
        atom.idx = i
    datafile.atoms = new_atoms

    datafile.write(args.out)
