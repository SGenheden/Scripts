# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to make a plummed-common input file for umbrella
sampling of the distance between a membrane and one or more solutes

The atom indices are taken from a Gromacs index file

Examples:
  make_plumed.py --solutes AAC1 AAC2
"""

import argparse

from sgenlib import groups

if __name__ == '__main__' :

    parser = argparse.ArgumentParser(description="Making plummed-common input")
    parser.add_argument('-n','--ndx',help="the input index file",default="index.ndx")
    parser.add_argument('-m','--mem',help="the group name of the membrame",default="Membrane")
    parser.add_argument('-s','--solutes',nargs="+",help="the group name of the solutes",default=[])
    args = parser.parse_args()

    ndxgroups = groups.read_ndxfile(args.ndx)

    for i, solute in enumerate(args.solutes, 1) :
        print "c%d: COM ATOMS=%d-%d"%(i, ndxgroups[solute][0], ndxgroups[solute][-1])
    print "c%d: COM ATOMS=%d-%d"%(len(args.solutes)+1, ndxgroups[args.mem][0],
                                    ndxgroups[args.mem][-1])
    print ""
    for i in range(1, len(args.solutes)+1):
        print "cv%d: DISTANCE ATOMS=c%d,c%d COMPONENTS"%(i, i, len(args.solutes)+1)
