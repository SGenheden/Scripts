# Author: Samuel Genheden, samuel.genheden@gmail.com

"""
Program to optimize LJ parameters

Used in a project to paramettrize ELBA ions
"""

import argparse

import numpy as np
import matplotlib.pylab as plt
import scipy.optimize as opt

from sgenlib import parsing

def objfunc(x,hfecoeff,iodcoeff,hfefit,iodfit):

    x2 = x*x
    x3 = x2*x
    hfex = hfecoeff[0]*x3+hfecoeff[1]*x2+hfecoeff[2]*x+hfecoeff[3]
    iodx = iodcoeff[0]*x3+iodcoeff[1]*x2+iodcoeff[2]*x+iodcoeff[3]
    hfediff = np.abs(np.divide(hfex - hfefit,hfefit)).sum()
    ioddiff = np.abs(np.divide(iodx - iodfit,iodfit)).sum()

    return hfediff + ioddiff

if __name__ == '__main__' :

    parser = argparse.ArgumentParser(description="Optimize LJ parameters")
    parser.add_argument('-s','--scanfile',help="the data from the scans")
    parser.add_argument('-f','--fitfile',help="the data points to fit")
    args = parser.parse_args()

    scans = parsing.parse2ndarray(args.scanfile)
    fitpnt = parsing.parse2ndarray(args.fitfile)

    hfecoeff = np.polyfit(scans[:,0],scans[:,1],3)
    iodcoeff = np.polyfit(scans[:,0],scans[:,2],3)

    xmin = scans[:,0].min()
    xmax = scans[:,0].max()
    x0 = np.ones(fitpnt.shape[0])*scans[:,0].mean()
    arglst = (hfecoeff,iodcoeff,fitpnt[:,0],fitpnt[:,1])
    bnds = np.ones((fitpnt.shape[0],2))
    bnds[:,0] *= xmin
    bnds[:,1] *= xmax

    res = opt.minimize(objfunc,x0,args=arglst,bounds=bnds)
    optx = res.x

    print "%7s\t%7s\t%7s\t%8s\t%8s\t%5s\t%5s"%("Rmin/2","Sigma","Eps","HFE(x)","HFE(f)","IOD(x)","IOD(f)")
    for x,fitvals in zip(optx,fitpnt):
        x2 = x*x
        x3 = x2*x
        hfex = hfecoeff[0]*x3+hfecoeff[1]*x2+hfecoeff[2]*x+hfecoeff[3]
        iodx = iodcoeff[0]*x3+iodcoeff[1]*x2+iodcoeff[2]*x+iodcoeff[3]
        print "%.5f\t%.5f\t%.5f\t%.3f\t%.3f\t%.3f\t%.3f"%(x,0.5**(1.0/6.0)*2*x,10**(-57.36*np.exp(-2.471*x)),fitvals[0],hfex,fitvals[1],iodx)
