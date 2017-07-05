


import argparse
import subprocess
import os
import sys
import tempfile

import numpy as np


if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description="Script to make tabulated potentials")
    argparser.add_argument('-e','--electrostatics',choices=["cut","switch","rf","zero"],help="the electrostatic force field",default="cut")
    argparser.add_argument('-v','--vanderwaals',choices=["cut","switch","zero"],help="the van der Waals force field",default="cut")
    argparser.add_argument('--cutoff',type=float,help="the cut-off",default=1.0)
    argparser.add_argument('--longoff',type=float,help="the long-range cut-off")
    argparser.add_argument('--switch',type=float,help="the switch distance",default=0.0)
    argparser.add_argument('--dielectric',type=float,help="the dielectric constant",default=1.0)
    argparser.add_argument('--epsilon',type=float,help="the dielectric constant",default=1.0)
    argparser.add_argument('-o', '--out', help="the output filename", default="table.xvg")
    args = argparser.parse_args()

    if args.longoff is None : args.longoff = args.cutoff

    cmd = os.path.join(os.path.dirname(sys.argv[0]), "make_tables ")

    tmpfile, tmpname = tempfile.mkstemp()
    with subprocess.Popen(cmd, shell=True,
        stdin=subprocess.PIPE, stdout=tmpfile, stderr=tmpfile).stdin as p :
        if args.electrostatics == "cut" :
            p.write("cut\n")
        elif args.electrostatics == "switch" :
            p.write("switched\n")
        elif args.electrostatics == "rf" :
            p.write("rf\n")
        elif args.electrostatics == "zero" :
            p.write("zero\n")

        if args.vanderwaals == "cut" :
            p.write("cut\n")
        elif args.vanderwaals == "switch" :
            p.write("switched\n")
        elif args.vanderwaals == "zero" :
            p.write("zero\n")

        p.write("%.2f\n"%args.cutoff)
        p.write("%.2f\n"%args.longoff)
        p.write("%.2f\n"%args.dielectric)

        if args.electrostatics == "rf" :
            p.write("%.2f\n"%args.epsilon)

        if args.electrostatics == "switch" or args.vanderwaals == "switch" :
            p.write("%.2f\n"%args.switch)

        p.write("%s\n"%args.out)

    os.remove(tmpname)
