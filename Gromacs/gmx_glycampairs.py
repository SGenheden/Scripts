# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to make proper Glycam pair parameters for 1-4 neighbours that
acpype.py cannot hanndle
"""

import argparse
import math

from sgenlib import gmx

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Program make Glycam pair params")
    parser.add_argument('-f','--file',help="the input top-file")
    parser.add_argument('-m','--molecule',help="the molecule name")
    args = parser.parse_args()

    reftop = gmx.TopFile(args.file)

    atypes = {}
    for atype in reftop.atomtypes :
        atypes[atype.name] = atype

    refmol = None
    for mol in reftop.moleculetypes :
        if mol.name == args.molecule :
            refmol = mol
            break
    if refmol is None :
        print "Could not find moleculetype %s. Exit."%args.molecule
        quit()

    for pair in refmol.pairs :
        atom1 = refmol.atoms[pair.atoms[0]-1]
        atom2 = refmol.atoms[pair.atoms[1]-1]

        ch1 = atom1.charge
        eps1 = atypes[atom1.type].epsilon
        sig1 = atypes[atom1.type].sigma
        ch2 = atom2.charge
        eps2 = atypes[atom2.type].epsilon
        sig2 = atypes[atom2.type].sigma
        # LB combination rule
        eps = math.sqrt(eps1*eps2)
        sig = 0.5 * (sig1 + sig2)
        
        print "%5d %5d  2   1.00   %15.6f %15.6f %15.6e %15.6e"%(
                pair.atoms[0], pair.atoms[1], ch1, ch2, sig, eps)
