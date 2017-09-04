# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to run external commands on a solutes from the Minnesota solvation database

Examples:
run_commands.py -db MNSol_alldata.txt -solvent hexanol -solutes hexanolwater.txt
                -commands insert_command
"""

import argparse

import dblib
from sgenlib import ambertools

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description="Script to run a command for all solutes in a list")
    argparser.add_argument('-db', '--db', help="the molecule database")
    argparser.add_argument('-solvent', '--solvent', help="the solvent")
    argparser.add_argument('-solutes','--solutes',help="the list of solutes")
    argparser.add_argument('-commands','--commands',help="a file with commands")
    args = argparser.parse_args()

    db = dblib.SolvDb(filename=args.db,type="abs",filehandle="^0")
    solutes = [s.strip() for s in open(args.solutes,'r').readlines()]
    commands = [s.strip() for s in open(args.commands,'r').readlines()]

    # Loop over all the database entries in the solute lists
    nerrors = 0
    for i,entry in enumerate(db.itersolutelist(args.solvent,solutes)):
        for cmd in commands :
            if i == 0 : print cmd.replace("$$",entry.FileHandle)
            try :
                ambertools.run_program("cmd",cmd.replace("$$",entry.FileHandle))
            except :
                print "!",
                nerrors += 1;
        if nerrors > 10 :
            raise Exception("Too many errors!")
        print entry.SoluteName, entry.FileHandle
