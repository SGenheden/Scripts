# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to fit a dihedral potential

The input file consist of rows with the following format
    DIHEDRAL MM QM
where DIHEDRAL is in degrees and MM and QM are the MM and QM energies, respectively

A three term dihedral potential will be fit to the difference between QM and MM (QM-MM)

Initial guess of shift1, k1 etc are to be given as the second argument

Uses in membrane engineering project

Examples
--------
fit_dihedral.py  mm-cc.dat "0 5 180 5 0 5"
"""
import sys

import numpy as np
from scipy.optimize import curve_fit

from sgenlib import parsing

data = parsing.parse2ndarray(sys.argv[1])

def dihedral_func(val, shift1, k1, shift2, k2, shift3, k3) :

    return k1*(1+np.cos(np.radians(1*x-shift1))) + \
            k2*(1+np.cos(np.radians(2*x-shift2))) + \
            k3*(1+np.cos(np.radians(3*x-shift3)))

x=data[:,0] # Dihedrals
y=data[:,2]-data[:,1] # QM - MM
popt, pcov = curve_fit(dihedral_func, x, y, map(float, sys.argv[2].split())) 

print "shift1 %7.2f"%popt[0]
print "k1 %10.3f"%popt[1]
print "shift2 %7.2f"%popt[2]
print "k2 %10.3f"%popt[3]
print "shift3 %7.2f"%popt[4]
print "k3 %10.3f"%popt[5]

for xi, yfit in zip(x, dihedral_func(x, popt[0], popt[1], popt[2], popt[3], popt[4], popt[5])) :
    print "%7.2f %10.3f"%(xi, yfit)
