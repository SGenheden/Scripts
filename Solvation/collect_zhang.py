# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to collect TI results for Zhang solutes

Examples:
collect_zhange.py -solvent octanol -solutes list
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
    argparser.add_argument('-solvent', '--solvent', help="the solvent")
    argparser.add_argument('-solutes','--solutes',help="the list of solutes")
    argparser.add_argument('-outdir','--outdir',help="the directory with output files",default=".")
    args = argparser.parse_args()

    solutes = [s.strip() for s in open(args.solutes,'r').readlines()]

    for solute in solutes:
        filename = os.path.join(args.outdir,
                    "out.dPotEngSS_%s_%s"%(solute,args.solvent))
        dglist = [_integrate(filename)]
        for repeat in range(2,5):
            filename = os.path.join(args.outdir,"R%d"%repeat,
                        "out.dPotEngSS_%s_%s"%(solute,args.solvent))
            try:
                dg_repeat = _integrate(filename)
            except:
                nomore = True
            else:
                dglist.append(dg_repeat)
                nomore = False
            if nomore: break
        if len(dglist) == 1:
            dg = dglist[0]
            std = 0.0
        else:
            dg = np.asarray(dglist).mean()
            std = np.asarray(dglist).std()/np.sqrt(len(dglist))
        print "%s\t%.3f\t%.3f"%(solute,-dg*4.184,std*4.184)
