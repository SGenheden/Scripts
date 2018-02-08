# Author: Samuel Genheden samuel.genheden@gmail.com
"""
Script to calculate the release free energy. The WRelease function is taken
from the offical APR code in september 2017.
"""

import sys

import numpy as np

def WRelease(kr,r0,kt,t0,kb,b0,T):
    """
       Do the analytical "release" of guest to standard concentration
    """
    ### SETUP
    # Example: print WRelease(5.0,24.00,100.0,180.0,100.0,180.0,298.15)
    # kr, r0: distance restraint (r) force constant, target value
    # kt, t0: angle restraint (theta) force constant, target value
    # kb, b0: angle restraint (b) force constant, target value
    # T: Temperature
    R = 1.987204118e-3 # kcal/mol-K, a.k.a. boltzman constant
    beta = 1/(T*R)
    rlb,rub,rst = [0.0,100.0,0.0001]  # r lower bound, upper bound, step size
    tlb,tub,tst = [0.0,np.pi,0.00005] # theta ",       ",            "
    blb,bub,bst = [0.0,np.pi,0.00005] # b     ",       ",            "
    def fr(val):
      return (val**2)*np.exp(-beta*kr*(val-r0)**2)
    def ft(val):
      return np.sin(val)*np.exp(-beta*kt*(val-np.radians(t0))**2)
    def fb(val):
      return np.sin(val)*np.exp(-beta*kb*(val-np.radians(b0))**2)
    ### Integrate
    rint,tint,bint = [0.0,0.0,0.0]
    intrange = np.arange(rlb,rub,rst)
    rint = np.trapz(fr(intrange),intrange)
    intrange = np.arange(tlb,tub,tst)
    tint = np.trapz(ft(intrange),intrange)
    intrange = np.arange(blb,bub,bst)
    bint = np.trapz(fb(intrange),intrange)
    return R*T*np.log(np.pi*(1.0/1660.0)*rint*tint*bint)

if __name__ == "__main__":


    lines = []
    with open(sys.argv[1], 'r') as f :
        lines = f.readlines()

    angfac = (180.0 / np.pi) * (180.0 / np.pi)
    for i, line in enumerate(lines) :
        if line.find("colvars d1") > -1 :
            kr = 0.5 * float(lines[i+1].strip().split()[1])
            r0 = float(lines[i+2].strip().split()[1])
        elif line.find("colvars a1") > -1 :
            kt = 0.5 * angfac * float(lines[i+1].strip().split()[1])
            t0 = float(lines[i+2].strip().split()[1])
        elif line.find("colvars a2") > -1 :
            kb = 0.5 * angfac * float(lines[i+1].strip().split()[1])
            b0 = float(lines[i+2].strip().split()[1])
    print "kr=%.3f r0=%.3f kt=%.3f t0=%.3f kb=%.3f b0=%.3f"%(kr,r0,kt,t0,kb,b0)
    print "Wrelease = %.3f"%WRelease(kr,r0,kt,t0,kb,b0,298)
