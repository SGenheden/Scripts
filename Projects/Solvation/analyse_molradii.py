# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to analyse the molecular radius of all solutes in a list
The solutes coordinates are taken from the Minnesota solvation database
"""

import argparse
import os

import numpy as np

import dblib

def _calc_radii(filename):

    with open(filename,"r") as f :
        data = [s.strip().split() for s in f.readlines()[3:]]
        data = np.array(data,dtype=float)[:,1:]
        len = data.max(axis=0)-data.min(axis=0)
        return np.max(len)
        

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description="Script to calculate molecular radii of all solutes in a list")
    argparser.add_argument('-db', '--db', help="the molecule database")
    argparser.add_argument('-solvent', '--solvent', help="the solvent")
    argparser.add_argument('-solutes','--solutes',help="the list of solutes")
    argparser.add_argument('-xyzdir','--xyzdir',help="the directory with all xyz files")
    args = argparser.parse_args()

    db = dblib.SolvDb(args.db,type="abs",filehandle="^0")
    solutes = [s.strip() for s in open(args.solutes,'r').readlines()]

    maxrad = -1000
    maxradsol = ""
    for entry in db.itersolutelist(args.solvent,solutes):
        rad = _calc_radii(os.path.join(args.xyzdir,entry.FileHandle+".xyz"))
        print "%-30s\t%.3f"%(entry.SoluteName,rad)
        if rad > maxrad :
            maxrad = rad
            maxradsol = entry.SoluteName

    print "\nMax:\n%-30s\t%.3f"%(maxradsol,maxrad)
    