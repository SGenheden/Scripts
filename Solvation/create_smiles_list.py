# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to list the SMILES string for each solute in a list
"""

import argparse
import os

from chemspipy import ChemSpider

import dblib

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description="Script to list SMILES for a solutes in a list")
    argparser.add_argument('-db', '--db', help="the molecule database")
    argparser.add_argument('-solvent1', '--solvent1', help="the first solvent", default="water")
    argparser.add_argument('-solvent2', '--solvent2', help="the second solvent", default="water")
    argparser.add_argument('-solutes','--solutes',help="the list of solutes")
    args = argparser.parse_args()

    db = dblib.SolvDb(filename=args.db,type="abs",filehandle="^0")
    solutes = [s.strip() for s in open(args.solutes,'r').readlines()]

    n = 1
    for entry1, entry2 in db.itersoluteoverlap(args.solvent1,args.solvent2):
        if  entry1.SoluteName in solutes and os.path.exists(entry1.FileHandle+".smi") :
            with open(entry1.FileHandle+".smi","r") as f :
                smi = f.readline().strip()
            print "%d %s %s %.3f"%(n, entry1.FileHandle, smi,
                    float(entry2.DeltaGsolv)-float(entry1.DeltaGsolv))
            n += 1
