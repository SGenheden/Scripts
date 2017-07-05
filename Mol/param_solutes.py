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
import parmed

from sgenlib import ambertools

leapcmd="""
source leaprc.gaff
loadamberprep %s.prepi
loadamberparams %s.frcmod
x=loadpdb %s.pdb
saveamberparm x %s_amb.top %s.crd
savepdb x %s_leap.pdb
quit
"""

def _make_resname(solute):
    if len(solute) == 3 :
        return solute
    elif len(solute) == 2 :
        solute=solute+"x"
    elif len(solute) == 1 :
        solute=solute+"xx"
    resname = ""
    pos = 0
    while len(resname) < 3:
        if solute[pos] in string.ascii_lowercase:
            resname += solute[pos]
        pos += 1
    return resname

def _write_gro(solute) :
    amb_top = parmed.amber.AmberParm("%s_amb.top"%solute, xyz="%s.crd"%solute)
    gmx_top = parmed.gromacs.GromacsTopologyFile.from_structure(amb_top)
    gmx_top.write("%s_gmx.top"%solute)

def _write_pdb(xyzname, resname) :

    bname = os.path.basename(xyzname)

    lines = []
    with open(xyzname, "r") as f :
        lines = f.readlines()[2:]

    s = parmed.Structure()
    r = parmed.Residue(name=resname, number=1, list=s.residues)
    s.residues.append(r)

    for i, atom_line in enumerate(lines, 1) :
        data = atom_line.strip().split()
        a = parmed.Atom(list=s.atoms,
                            atomic_number=parmed.periodic_table.AtomicNum[data[0]],
                            name=data[0]+"%d"%i)
        a.number = i
        a.xx = float(data[1])
        a.xy = float(data[2])
        a.xz = float(data[3])
        r.add_atom(a)
        s.atoms.append(a)

    s.write_pdb(os.path.splitext(bname)[0]+".pdb")

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description="Script to parametrize all solutes in a list")
    argparser.add_argument('filename',help="the filename of a list")
    argparser.add_argument('-xyzdir','--xyzdir',help="the directory with all xyz files",default=".")
    argparser.add_argument('--noclean',action="store_true",help="remove unnecessary files",default=False)
    argparser.add_argument('--postfix',help="a postfix for the xyz filenames",default="")
    args = argparser.parse_args()

    solutes = [s.strip() for s in open(args.filename,'r').readlines()]

    n = 0
    for solute0 in solutes:
        solute = solute0+args.postfix
        # If the prepi file exists, we assume this solute has been done already
        if not os.path.exists(solute+".prepi"):
            resname = _make_resname(solute0.lower())
            # if not, write out a pdb, make an Amber prepi and frcmod file
            _write_pdb(os.path.join(args.xyzdir,solute+".xyz"),resname)
            hasparm = True
            try :
                ambertools.run_antechamber(solute+".pdb",0,resnam=resname)
            except :
                print "! ",
                hasparm = False
            if hasparm :
                ambertools.run_parmchk(solute+".pdb")
                # then make an Amber prmtop and prmcrd file
                with open("leapcom","w") as f :
                    f.write(leapcmd%tuple([solute]*6))
                ambertools.run_program("tleap","tleap -f leapcom")
                _write_gro(solute)
            if not args.noclean  :
                os.remove(solute+".pdb") # Remove unnecessary pdb file
            if not os.path.exists(solute+".prepi"): print "**"

        print solute0
        n += 1

    print "Looped over %d solutes"%n
