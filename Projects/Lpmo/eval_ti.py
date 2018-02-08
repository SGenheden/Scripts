# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to collect TI output
"""

import argparse

import numpy as np

from sgenlib import parsing

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description="Script to collect output and do TI")
    argparser.add_argument('outlist',nargs="+", help="the output list")
    argparser.add_argument('-dc','--deriv-col',type=int)
    args = argparser.parse_args()

    skip = 1.0 / 3.0

    dvdl = np.zeros((21,4))
    dvdl[:, 0] = np.arange(0.0,1.05,0.05)
    for i, filename in enumerate(args.outlist) :
        data = parsing.parse2ndarray(filename)
        start = int(len(data)*skip)
        dvdl[i,1] = data[start:, 1].mean()
        dvdl[i,2] = data[start:, 2].mean()
        dvdl[i,3] = data[start:, 3].mean()
    dvdl[-1,3] =  np.polyfit(dvdl[-3:-1,0],dvdl[-3:-1,3],1).sum()
    dg = np.trapz(dvdl[:,1],x=dvdl[:,0]) + np.trapz(dvdl[:,2],x=dvdl[:,0]) + np.trapz(dvdl[:,3],x=dvdl[:,0])
    print "dG = %.3f"%dg
