# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to analyse the overlap of solutes in two solvents
The solutes are taken from the Minnesota solvation database
"""

import argparse
import re

import dblib

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description="Script to analyse overlap of solutes in two solvents")
    argparser.add_argument('-db', '--db', help="the molecule database")
    argparser.add_argument('-solvent1', '--solvent1', help="the first solvent")
    argparser.add_argument('-solvent2', '--solvent2', help="the second solvent")
    argparser.add_argument('-list','--list', help="list of solutes")
    args = argparser.parse_args()

    db = dblib.SolvDb(filename=args.db,type="abs",filehandle="^0")

    list0 = []
    with open(args.list, "r") as f :
        list0 = [l.strip() for l in f.readlines()]

    dblist = []
    for entry1, entry2 in db.itersoluteoverlap(args.solvent1,args.solvent2) :
        if entry1.SoluteName not in list0 :
            print entry1.SoluteName
        dblist.append(entry1.SoluteName)

    print "---"

    with open("newlist.txt", "w") as f :
        for l in list0 :
            if l not in dblist :
                print l
            else :
                f.write("%s\n"%l)
