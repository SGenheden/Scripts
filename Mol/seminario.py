# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to calculate bond force constant with the Seminario method
for specific bonds

Ref:
JM Seminario, Int. J. Quant. Chem. 1996, 60(7), 1271-1277

Examples:
seminario.py -f model.fchk -s model.gro -b 1:N-4:CU 1:ND1-4:CU 2:NE2-4:CU
"""

import argparse
import os

import numpy as np

from sgenlib import pdb

def _seminario(atm1, atm2, hess, scaling) :
    """
    Calculate bond force constant by analysing the Hessian submatrix
    Average over atm1-atm2 and atm2-atm1 force constants
    """
    def _calc_force(atm1, atm2, hess, vec12):
        submat = -hess[3*atm1.idx:3*atm1.idx+3, 3*atm2.idx:3*atm2.idx+3]
        eigval, eigvec = np.linalg.eig(submat)
        fc = 0.0
        for i in range(3):
            fc += eigval[i] * np.abs(np.dot(eigvec[:,i], vec12))
        return fc

    # Vector from atm1 to atm2
    vec12 = atm1.xyz - atm2.xyz
    vec12 /= np.sqrt(atm1.distance2(atm2))

    f12 = _calc_force(atm1, atm2, hess, vec12)
    f21 = _calc_force(atm2, atm1, hess, vec12)
    # 2240.87 is from Hartree/Bohr ^2 to kcal/mol/A^2
    # 418.4 is kcal/mol/A^2 to kJ/mol/nm^2
    return scaling * 2240.87 * 418.4 * 0.5 * (f12+f21)

def _find_pair(struct, bond):
    """
    Parse atom specification to Atom objects
    """
    atm1, atm2 = bond.split("-")
    res1, atm1 = atm1.split(":")
    res2, atm2 = atm2.split(":")

    a1 = struct.residues[int(res1)-1].atom_by_name(atm1)
    a2 = struct.residues[int(res2)-1].atom_by_name(atm2)
    return a1, a2

def _parse_fchk(filename):
    """
    Parse Gaussian09 formated checkpoint file for coordinates
    and Hessian
    """
    def _parse_array(f, startline, endline):
        arr = []
        line = "None"
        while line and not line.startswith(startline) :
            line = f.readline()
        while line and not line.startswith(endline):
            line = f.readline()
            if not line.startswith(endline):
                arr.extend(line.strip().split())
        return np.array(arr, dtype=float)

    crds = None
    hess = None
    with open(filename, "r") as f :
        # First the coordinates
        crds = _parse_array(f, 'Current cartesian coordinates', 'Force Field')
        # Then the Hessian
        hess = _parse_array(f, 'Cartesian Force Constants', 'Dipole Moment')

    # Make the Hessian in square form
    i = 0
    n = len(crds)
    hess_sqr = np.zeros([n, n])
    for j in range(n):
        for k in range(j+1):
            hess_sqr[j,k] = hess[i]
            hess_sqr[k,j] = hess[i]
            i += 1

    return crds, hess_sqr

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description="Script to compute bond force constant with Seminario")
    argparser.add_argument('-f', '--checkpoint', help="the format checkpoint file")
    argparser.add_argument('-s', '--struct', help="a gro or pdb file with atom/residue information")
    argparser.add_argument('-b', '--bonds', nargs="+", help="the bonds")
    argparser.add_argument('--scaling', type=float, help="the frequency scaling factor", default=0.963)
    args = argparser.parse_args()

    crds, hess = _parse_fchk(args.checkpoint)

    struct_orig = pdb.PDBFile(args.struct)
    struct = pdb.PDBFile(args.struct)
    for i, atom in enumerate(struct.atoms, 1):
        # 0.529177249 is Bohr to A
        atom.set_xyz(crds[3*i-3:3*i]*0.529177249)
    base, ext = os.path.splitext(args.struct)
    struct.write(base+"_opt"+ext)

    print "Bond\tr(x-ray)\tr(opt) [nm]\tk [kJ/mol/nm2]"
    for bond in args.bonds :
        atm1, atm2 = _find_pair(struct, bond)
        dist = np.sqrt(atm1.distance2(atm2))

        atm1_orig = struct_orig.atoms[atm1.idx]
        atm2_orig = struct_orig.atoms[atm2.idx]
        dist_orig = np.sqrt(atm1_orig.distance2(atm2_orig))

        force = _seminario(atm1, atm2, hess, args.scaling)

        print "%s\t%d\t%d\t%.4f\t%.4f\t%.4f"%(bond, atm1.idx+1, atm2.idx+1, dist_orig*0.1, dist*0.1, force)
