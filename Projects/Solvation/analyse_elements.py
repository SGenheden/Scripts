# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to analyse the elements of all solutes in a list
The solutes are taken from the Minnesota solvation database
"""

import argparse
import re

import dblib

def _elements(form) :
    """
    Parse elements from molecular formula
    """
    elements = []
    for part in re.findall("[A-Z]+[0-9]+",form):
        m = re.match("([A-Z]+)([0-9]+)",part)
        elements.append(m.group(1))
    return set(elements)

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description="Script to analyse the elements of solutes")
    argparser.add_argument('-db', '--db', help="the molecule database")
    argparser.add_argument('-solvent', '--solvent', help="the solvent")
    argparser.add_argument('-solutes','--solutes',help="the list of solutes")
    args = argparser.parse_args()

    db = dblib.SolvDb(filename=args.db,type="abs",filehandle="^0")
    solutes = [s.strip() for s in open(args.solutes,'r').readlines()]
    gaffelements = "C H F CL BR I N O P S".split()

   # Loop over all the database entries in the solute lists
    for entry in db.itersolutelist(args.solvent,solutes):
        if not _elements(entry.Formula).issubset(gaffelements) :
            print entry.SoluteName,entry.Formula,entry.FileHandle