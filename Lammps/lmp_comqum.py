# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to setup a ComQum-ELBA run
Part of ComQum-ELBA
"""

import argparse

import numpy as np

from sgenlib import lammps

control_str = \
"""$point_charges  file=pointcharges
$operating system unix
$path
$symmetry c1
$coord   file=coord
$scfintunit
 unit=30       size=0        file=twoint
$scfconv   6
$scfiterlimit   300
$mvd
$moments
$pop
$statistics off
$optimize
   internal   off
   cartesian  on
$forceapprox   file=forceapprox
$interconversion  off
  qconv=1.d-10
   maxiter=25
$coordinateupdate
   dqmax=0.3
   interpolate  on
   statistics    5
$forceinit on
   diag=default
$forceupdate
   ahlrichs numgeo=0 mingeo=3 maxgeo=4 mode=<g|dq> dynamic fail=0.3
   threig=0.005  reseig=0.005 thrbig=3.0 scale=1.00 damping=0.0
$grad   file=gradient
$suspend off
$lock off
$dft
   functional b3-lyp
   gridsize   m3
$end
"""

mass2element = {1.00800 : "h", 12.01000 : "c", 16.00000 : "o"}
a2au = 1.0 / 0.529177249

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Prepare a ComQum-ELBA run from a LAMMPS datafile")
    parser.add_argument('file', help="the data file.")
    parser.add_argument('-m', '--mol', type=int, help="the molecule id of the QM system")
    args = parser.parse_args()

    if args.file is None:
        print "No input file specified. Exiting!"
        quit()

    datafile = lammps.Datafile(args.file)
    qmatoms = [atom for atom in datafile.atoms if atom.molecule == args.mol]
    other_atoms = [atom for atom in datafile.atoms if atom.molecule != args.mol]
    atom_ids = [atom.idx for atom in qmatoms]

    with open("comqum.dat", "w") as f :
        f.write("$title\nNew ComQum-ELBA project\n")
        f.write("$protein fixed\n")
        f.write("$correspondance_list\n")
        f.write("%6d\n"%len(qmatoms))
        for i, id in enumerate(atom_ids, 1) :
            f.write("%6d "%id)
            if i % 12 == 0 :
                f.write("\n")
        f.write("\n")
        f.write("$end\n")

    with open("control", "w") as f :
        f.write(control_str)

    with open("coord", "w") as f :
        f.write("$coord\n")
        for atom in qmatoms :
            f.write("%16.11f %16.11f %16.11f %s\n"%(atom.x*a2au, atom.y*a2au, atom.z*a2au, mass2element[atom.density]))
        f.write("$end\n")

    with open("pointcharges", "w") as f :
        f.write("$point_charges mxrank=1\n")
        for atom in other_atoms :
            f.write("%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n"%(
                atom.x*a2au, atom.y*a2au, atom.z*a2au,
                atom.q,
                atom.mux*a2au, atom.muy*a2au, atom.muz*a2au,
            ))
        f.write("$end\n")
