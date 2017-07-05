# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to convert a Gromacs top file to an itp-file.

This is not a bullet-proof procedure that works in all situations!
Works for small solutes created with param_solutes.py and this was the usage
it was created for

Examples:
top2itp.py file1.itp file2.itp
"""

import argparse
import os

def _convert(filename, atomtypes) :

    lines = []
    with open(filename, "r") as f :
        lines = f.readlines()

    with open(os.path.splitext(filename)[0]+".itp", "w") as f :
        i = 0
        while i < len(lines) and lines[i].find("defaults") == -1 :
            f.write(lines[i])
            i += 1
        while i < len(lines) and lines[i].find("atomtypes") == -1 :
            i += 1
        i += 2
        while i < len(lines) and len(lines[i].strip()) > 5 :
            atype = lines[i].strip().split()[0]
            atomtypes[atype] = lines[i]
            i += 1
        while i < len(lines) and lines[i].find("system") == -1 :
            f.write(lines[i])
            i += 1

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description="Script convert a Gromacs top file to an itp-file")
    argparser.add_argument('filenames',nargs="+",help="the filename of a list")
    args = argparser.parse_args()

    atomtypes = {}
    for filename in args.filenames :
        print "Converting %s..."%filename
        _convert(filename, atomtypes)

    print "These are the atom types that need to be inserted elsewhere:"
    for at in sorted(atomtypes.keys()) :
        print atomtypes[at],
