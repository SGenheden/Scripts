# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to parametrise solutes from a list

Examples:
param_solutes.py list
"""

import argparse
import os
import string

import numpy as np
import MDAnalysis as md

from sgenlib import ambertools

leapcmd="""
source leaprc.gaff
loadamberprep %s.prepi
loadamberparams %s.frcmod
x=loadpdb %s.pdb
saveamberparm x %s.top %s.crd
savepdb x %s_leap.pdb
quit
"""

def _make_resname(solute):
    if len(solute) == 3 : return solute
    resname = ""
    pos = 0
    while len(resname) < 3:
        if solute[pos] in string.ascii_lowercase:
            resname += solute[pos]
        pos += 1
    return resname

def _write_pdb(xyzname, resname) :

    bname = os.path.basename(xyzname)

    # Read it in to a MDAnalysis universe
    reader = md.Universe(bname,bname)
    for i,atom in enumerate(reader.selectAtoms("all"),1) :
        atom.resname = resname
        atom.name = atom.name+"%d"%i

    # Write it out as a pdb
    writer = md.coordinates.PDB.PDBWriter(os.path.splitext(bname)[0]+".pdb", universe=reader)
    writer.write(reader.trajectory.ts)
    writer.close()

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description="Script to parametrize all solutes in a list")
    argparser.add_argument('filename',help="the filename of a list")
    argparser.add_argument('-xyzdir','--xyzdir',help="the directory with all xyz files",default=".")
    args = argparser.parse_args()

    solutes = [s.strip() for s in open(args.filename,'r').readlines()]

    n = 0
    for solute in solutes:
        # If the prepi file exists, we assume this solute has been done already
        if not os.path.exists(solute+".prepi"):
            resname = _make_resname(solute.lower())
            # if not, write out a pdb, make an Amber prepi and frcmod file
            _write_pdb(os.path.join(args.xyzdir,solute+".xyz"),resname)
            ambertools.run_antechamber(solute+".pdb",0,resnam=resname)
            ambertools.run_parmchk(solute+".pdb")
            # then make an Amber prmtop and prmcrd file
            with open("leapcom","w") as f :
                f.write(leapcmd%tuple([solute]*6))
            ambertools.run_program("tleap","tleap -f leapcom")
            os.remove(solute+".pdb") # Remove unnecessary pdb file
            if not os.path.exists(solute+".prepi"): print "**"

        print solute
        n += 1

    print "Looped over %d solutes"%n
