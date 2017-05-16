# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to modify LJ parameters in LAMMPS inclusion file
"""

import argparse
from collections import namedtuple

import numpy as np

from sgenlib import lammps

if __name__ == '__main__':

    # Command-line input
    parser = argparse.ArgumentParser(description="Modifying LJ parameters in LAMMPS inclusion file")
    parser.add_argument('file', help="the inclusion file.")
    parser.add_argument('-t', '--type', type=int, help="the atom type to modify")
    parser.add_argument('-p', '--parameters', nargs="+",help="the new parameters")
    args = parser.parse_args()

    if args.file is None:
        print "No input file specified. Exiting!"
        quit()

    infile = lammps.Includefile(filename=args.file)

    Ljparam = namedtuple("Ljparam",['epsilon', 'sigma'])
    params = []
    if len(args.parameters) == 2:
        params.append(Ljparam(float(args.parameters[0]),float(args.parameters[1])))
    else:
        with open(args.parameters[0],"r") as f:
            for line in f.readlines():
                val = map(float,line.strip().split())
                params.append(Ljparam(*val))

    # Find all mixed pairs
    mixed_pairs = []
    for pair in infile.pair_coeff:
        if (pair.comment.find("AA-CG mixed") > -1 or pair.comment.find("Mixed using Lorentz-Berthelot rules") > -1):
            mixed_pairs.append(pair)
    lj_same = infile.lj_same_dict()

    for pi,param in enumerate(params,1):
        lj_same[args.type].epsilon = param.epsilon
        lj_same[args.type].sigma = param.sigma
        for pair in mixed_pairs:
            if args.type in [pair.jatom,pair.iatom] :
                e = lj_same[pair.jatom].epsilon if pair.iatom == args.type else lj_same[pair.iatom].epsilon
                s = lj_same[pair.jatom].sigma if pair.iatom == args.type else lj_same[pair.iatom].sigma
                pair.epsilon = np.sqrt(e*param.epsilon)
                pair.sigma = (s + param.sigma)/2.0

        # Write output
        infile.write(args.file+"_param%d"%pi,ljres=10)
