# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Script to calculate PMF properties, i.e. penetration and water/lipid barriers
"""

import sys

import numpy as np
from scipy.interpolate import interp1d

from sgenlib import parsing
from sgenlib import snum

data = parsing.parse2ndarray(sys.argv[1])
x = data[:,0]
pmf = data[:,1]
err = data[:,2]

mini = np.argmin(pmf)

dGdepth = pmf[mini] - pmf[-1]
dGpen   = pmf[0] - pmf[mini]

errdepth = np.sqrt(err[mini]**2 + err[-1]**2)
errpen   = np.sqrt(err[0]**2 + err[mini]**2)

kT = 300*0.00831446210
expav = np.exp(-pmf/kT)
expstd = np.abs(expav*(-err/kT))
bndint, bndstd = snum.trapz(x, expav, expstd)
freeint = np.trapz(np.ones(x.shape[0]), x)
dGbind = -kT * np.log(bndint/freeint)
errbind = np.abs(kT*bndstd/bndint)

if len(sys.argv) > 2 :
    water_dens = parsing.parse2ndarray(sys.argv[2])
    f = interp1d(water_dens[:,0],water_dens[:,1],kind="cubic")
    rho = f(x)
    rho = rho/rho.max()
    APL = 0.62
    const = APL*np.power(10,-24.0)/(786.14*1.66*np.power(10,-27.0))
    expav = np.exp(-pmf/kT) - rho
    expstd = np.abs(expav*(-err/kT))
    bnd,bndstd = snum.trapz(x*0.1, expav, expstd)
    if bnd < 0 :
      expstd[expstd<0] = 0
      expav[expav<0] = 0
      bnd,bndstd = snum.trapz(x, expav, expstd)
    part = np.log10(const*bnd)
    partstd = np.abs(bndstd*const/(bnd*np.log(10)))
    print "%.3f\t%.3f\t"%(part,partstd),



print "\t".join("%.3f"%g for g in (dGdepth, errdepth,
                                    dGpen, errpen, pmf[0], err[0], dGbind, errbind))
