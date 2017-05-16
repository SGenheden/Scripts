# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Strip water and displace bilayer
"""

import argparse

import numpy as np

from sgenlib import pdb

def _calc_dimensions(struct, sel, verbose) :

    sorti = np.argsort(struct.xyz[sel, 2])
    midi = int(0.5 * len(sel))
    low_len = np.floor(struct.xyz[sel, 2][sorti[:midi]].mean())
    upp_len = np.ceil(struct.box[2] - struct.xyz[sel, 2][sorti[midi:]].mean())
    watlen = np.ceil(0.5*(low_len + upp_len))
    if verbose :
        print "Low length = %d"%low_len
        print "Upp length = %d"%upp_len
        print "Watbox length = %d"%watlen
    return watlen

if __name__ == '__main__' :

    parser = argparse.ArgumentParser(description="Make a clean bilayer")
    parser.add_argument('-f','--file',help="the input file")
    parser.add_argument('-o','--out',help="the output file",default="bilayer.gro")
    parser.add_argument('-p','--patom',help="the name of the phosphor atom", default="P")
    parser.add_argument('-w','--watom',help="the name of the water atom", default="OH2")
    parser.add_argument('-s','--scale',type=float,help="the displacement scale", default=0.0)
    parser.add_argument('-v','--verbose',action="store_true", help="verbose output", default=False)
    args = parser.parse_args()

    struct = pdb.PDBFile(filename=args.file)

    residues = []
    atoms = []
    xyz = []
    for residue in struct.residues :
        if residue.atoms[0].name.strip() != args.watom :
            residues.append(residue)
            atoms.extend(residue.atoms)
            xyz.extend([atom.xyz for atom in residue.atoms])
    xyz = np.asarray(xyz)
    struct.xyz = xyz
    struct.atoms = atoms
    struct.residues = residues

    sel = np.asarray([i for i, atom in enumerate(struct.atoms) if atom.name.strip() == args.patom])
    watlen = _calc_dimensions(struct, sel, args.verbose)

    if args.scale > 0 :
        displ = args.scale * watlen
        for i, atom in enumerate(struct.atoms) :
            atom.z += displ
            atom.xyz[2] += displ
            xyz[i, 2] += displ
        struct.box[2] += 2*displ

        if args.verbose :
            print ""
        _calc_dimensions(struct, sel, args.verbose)

    sorti = np.argsort(struct.xyz[sel, 2])
    midi = int(0.5 * len(sel))
    if args.verbose :
        print " -- Bounding boxes -- "
    print "%.3f %.3f %.3f"%(0.0, 0.0, 0.0)
    print "%.3f %.3f %.3f"%(struct.box[0], struct.box[1],
                                struct.xyz[sel, 2][sorti[:midi]].mean())
    if args.verbose :
        print " -- --"
    print "%.3f %.3f %.3f"%(0.0, 0.0, struct.xyz[sel, 2][sorti[midi:]].mean())
    print "%.3f %.3f %.3f"%(struct.box[0], struct.box[1], struct.box[2])

    struct.write(args.out)
