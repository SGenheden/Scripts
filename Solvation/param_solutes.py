# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to parametrise solutes from the Minnesota solvation database

Examples:
param_solutes.py -db MNSol_alldata.txt -solvent hexanol -solutes hexanolwater.txt
                -xyzdir MNSolDatabase-v2012/all_solutes/
"""

import argparse
import os

import numpy as np
import MDAnalysis as md

import dblib
from sgenlib import ambertools

Element = [ 'EP',
            'H' ,'He','Li','Be','B' ,'C' ,'N' ,'O' ,'F' ,'Ne','Na','Mg',
            'Al','Si','P' ,'S' ,'Cl','Ar','K' ,'Ca','Sc','Ti','V' ,'Cr',
            'Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
            'Rb','Sr','Y' ,'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
            'In','Sn','Sb','Te','I' ,'Xe','Cs','Ba','La','Ce','Pr','Nd',
            'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf',
            'Ta','W' ,'Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po',
            'At','Rn','Fr','Ra','Ac','Th','Pa','U' ,'Np','Pu','Am','Cm',
            'Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs',
            'Mt','Ds','Rg','Cn','Uut','Uuq','Uup','Uuh','Uus','Uuo' ]

leapcmd="""
source leaprc.gaff
loadamberprep %s.prepi
loadamberparams %s.frcmod
x=loadpdb %s.pdb
saveamberparm x %s.top %s.crd
savepdb x %s_leap.pdb
quit
"""

def _write_pdb(xyzname) :

    bname = os.path.basename(xyzname)

    # Read in the databse xyz file
    xyz = None
    with open(xyzname,"r") as f :
        lines = f.readlines()
        data = [s.strip().split() for s in lines[3:]]
        data = np.array(data,dtype=float)

    # Write it out in the proper format
    with open(bname,"w") as f:
        f.write("%d\n"%data.shape[0])
        f.write("0\n")
        for pos in data :
            f.write("%s %15.8f %15.8f %15.8f\n"%(Element[int(pos[0])],pos[1],pos[2],pos[3]))

    # Read it in to a MDAnalysis universe
    resnam = os.path.splitext(bname)[0][-3:]
    reader = md.Universe(bname,bname)
    for i,atom in enumerate(reader.selectAtoms("all"),1) :
        atom.resname = resnam
        atom.name = atom.name+"%d"%i

    # Write it out as a pdb
    writer = md.coordinates.PDB.PDBWriter(os.path.splitext(bname)[0]+".pdb", universe=reader)
    writer.write(reader.trajectory.ts)
    writer.close()

    # Remove the temporary xyz file
    os.remove(bname)

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description="Script to parametrize all solutes in a list")
    argparser.add_argument('-db', '--db', help="the molecule database")
    argparser.add_argument('-solvent', '--solvent', help="the solvent")
    argparser.add_argument('-solutes','--solutes',help="the list of solutes")
    argparser.add_argument('-xyzdir','--xyzdir',help="the directory with all xyz files")
    args = argparser.parse_args()

    db = dblib.SolvDb(filename=args.db,type="abs",filehandle="^0")
    solutes = [s.strip() for s in open(args.solutes,'r').readlines()]

    # Loop over all the database entries in the solute lists
    n = 0
    for entry in db.itersolutelist(args.solvent,solutes):
        # If the PDB file exists, we assume this solute has been done already
        if not os.path.exists(entry.FileHandle+".prepi"):
            # if not, write out a pdb, make an Amber prepi and frcmod file
            _write_pdb(os.path.join(args.xyzdir,entry.FileHandle+".xyz"))
            ambertools.run_antechamber(entry.FileHandle+".pdb",0,resnam=entry.FileHandle[-3:])
            ambertools.run_parmchk(entry.FileHandle+".pdb")
            # then make an Amber prmtop and prmcrd file
            with open("leapcom","w") as f :
                f.write(leapcmd%tuple([entry.FileHandle]*6))
            ambertools.run_program("tleap","tleap -f leapcom")
            os.remove(entry.FileHandle+".pdb") # Remove unnecessary pdb file
            if not os.path.exists(entry.FileHandle+".prepi"): print "**"

        print entry.SoluteName
        n += 1

    print "Looped over %d solutes"%n
