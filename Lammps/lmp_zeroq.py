# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to zero charge of a QM system
Part of ComQum-ELBA
"""

import argparse

from sgenlib import lammps

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Truncate a LAMMPS datafile")
    parser.add_argument('file', help="the data file.")
    parser.add_argument('-m', '--mol', type=int, help="the molecule id of the QM system")
    parser.add_argument('-o', '--out',help="the output file")
    args = parser.parse_args()

    if args.file is None:
        print "No input file specified. Exiting!"
        quit()

    datafile = lammps.Datafile(args.file)
    for atom in datafile.atoms :
        if atom.molecule == args.mol :
            atom.q = 0.0

    if args.out is None :
        datafile.write(args.file)
    else :
        datafile.write(args.out)
