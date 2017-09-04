# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to obtain SMILES for the solutes in the Minnesota solvation database

Examples:
param_solutes.py -db MNSol_alldata.txt -solvent hexanol -solutes hexanolwater.txt
                -xyzdir MNSolDatabase-v2012/all_solutes/
"""

import argparse
import os

from chemspipy import ChemSpider

import dblib

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description="Script to obtain SMILES for a solutes in a list")
    argparser.add_argument('-db', '--db', help="the molecule database")
    argparser.add_argument('-solvent', '--solvent', help="the solvent", default="water")
    argparser.add_argument('-solutes','--solutes',help="the list of solutes")
    args = argparser.parse_args()

    db = dblib.SolvDb(filename=args.db,type="abs",filehandle="^0")
    solutes = [s.strip() for s in open(args.solutes,'r').readlines()]

    if os.getenv("SPIDERKEY") is None :
        print "SPIDERKEY environmental variable not set! Exit."
        quit()
    cs = ChemSpider(os.getenv("SPIDERKEY"))

    # Loop over all the database entries in the solute lists
    n = 0
    for entry in db.itersolutelist(args.solvent,solutes):
        if  os.path.exists(entry.FileHandle+".smi") : continue
        hits = cs.search(entry.SoluteName)
        if len(hits) > 0 :
            smi = hits[0].smiles
            with open(entry.FileHandle+".smi","w") as f :
                f.write("%s\n"%smi)
        else :
            print entry.SoluteName, entry.FileHandle
        n += 1

    print "Looped over %d solutes"%n
