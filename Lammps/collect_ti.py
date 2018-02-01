# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to collect TI results for a list of output

Examples:
collect_ti.py list
"""

import argparse
import os

import numpy as np

def _integrate(filename):
    """
    Read a PMF from an Lammps output file and integrate it
    """
    lines = []
    with open(filename,"r") as f :
        lines = f.readlines()
    lines.append("0.0 1.0 0.0")
    data = np.array([line.strip().split() for line in lines[2:]],dtype=float)
    lam = data[:,1]
    grad = data[:,2]
    # Linear interpolation to lambda=1 from the last two lambda values simulated
    grad[-1] = np.polyfit(lam[-3:-1],grad[-3:-1],1).sum()
    return np.trapz(grad,x=lam)

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description="Script to collect output and do TI")
    argparser.add_argument('outlist', help="the output list")
    argparser.add_argument('-outdir','--outdir',help="the directory with output files",default=".")
    argparser.add_argument('--fulloutput',action="store_true",default=False,help="if the output the individual results")
    argparser.add_argument('--prefix',help="the filename prefix")
    argparser.add_argument('--postfix',help="the filename postfix")
    args = argparser.parse_args()

    if args.prefix is not None:
        args.prefix = args.prefix+"-"
    else:
        args.prefix = ""

    if args.postfix is not None:
        args.postfix = "-"+args.postfix
    else:
        args.postfix = ""

    filenames = [s.strip() for s in open(args.outlist,'r').readlines()]

    for filename0 in filenames:
        filename = os.path.join(args.outdir,
                    "out.dPotEngSS_%s%s%s"%(args.prefix,filename0,args.postfix))
        dglist = [_integrate(filename)]
        for repeat in range(2,11):
            filename = os.path.join(args.outdir,"R%d"%repeat,
                        "out.dPotEngSS_%s%s%s"%(args.prefix,filename0,args.postfix))
            try:
                dg_repeat = _integrate(filename)
            except:
                nomore = True
            else:
                dglist.append(dg_repeat)
                nomore = False
            if nomore: break
        print "%s"%(filename0),
        if args.fulloutput :
            print "\t"+"\t".join("%.3f"%(-dg*4.184) for dg in dglist),
        if len(dglist) == 1:
            dg = dglist[0]
            std = 0.0
        else:
            dg = np.asarray(dglist).mean()
            std = np.asarray(dglist).std()/np.sqrt(len(dglist))
        print "\t%.3f\t%.3f"%(-dg*4.184,std*4.184)
