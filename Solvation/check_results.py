# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to
"""

import argparse
import os
import subprocess

import numpy as np

import dblib

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description="Script to check output")
    argparser.add_argument('-db', '--db', help="the molecule database")
    argparser.add_argument('-solvent', '--solvent', help="the solvent")
    argparser.add_argument('-solutes','--solutes',help="the list of solutes")
    argparser.add_argument('-outdir','--outdir',help="the directory with output files",default=".")
    args = argparser.parse_args()

    db = dblib.SolvDb(filename=args.db,type="abs",filehandle="^0")
    solutes = [s.strip() for s in open(args.solutes,'r').readlines()]

    # Loop over all the database entries in the solute lists
    reruns = []
    for i,entry in enumerate(db.itersolutelist(args.solvent,solutes),2):
        filename = os.path.join(args.outdir,"out.dPotEngSS_%s_%s"%(args.solvent,entry.FileHandle))
        try:
            with open(filename,"r") as f :
                lines = f.readlines()
            if float(lines[-1].strip().split()[1]) != 0.96:
                print "Simulation did not complete for %s = %s"%(entry.SoluteName,entry.FileHandle)
                reruns.append(entry.FileHandle)
        except :
            print "Unable to find output for %s = %s"%(entry.SoluteName,entry.FileHandle)
            subprocess.call("tail -n 1 screen.ti_%s_%s"%(args.solvent,entry.FileHandle), shell=True)
    if reruns :
        print "for X in %s"%" ".join(reruns)
