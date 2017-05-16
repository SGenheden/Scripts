# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to insert a solute at different depths in a liposome

Two solutes will be inserted by moving them in the opposite direction
from the center of the liposome

The solute will all likely overlap with the liposome and therefore
it is essential to grow it before attemping e.g. umbrella sampling

Example:
  insert_in_membrane.py -f hal.gro -b liposome.gro -r 120
"""

import argparse

import numpy as np
import numpy.random as random

from sgenlib import pdb
from sgenlib import geo

if __name__ == '__main__' :

    parser = argparse.ArgumentParser(description="Insert a solute in a liposome for umbrella sampling")
    parser.add_argument('-f','--file',help="the solute file")
    parser.add_argument('-b','--box',help="the liposome box")
    parser.add_argument('-o','--out',help="the output file",default="inserted")
    parser.add_argument('-r','--radii', nargs="+", type=float,help="the insertion radius from the liposome center")
    parser.add_argument('-p','--patom',help="the name of the phosphate atom", default="PO4")
    parser.add_argument('-n','--nvec',type=int,help="the index of a previously employed vector")
    args = parser.parse_args()

    rads_rev = args.radii[::-1]
    rads_pack = zip(args.radii,rads_rev)

    boxfile = pdb.PDBFile(filename=args.box)

    solutefile1 = pdb.PDBFile(filename=args.file)
    solutefile2 = pdb.PDBFile(filename=args.file)
    solxyz = np.array(solutefile1.xyz, copy=True)
    solcent = solxyz.mean(axis=0)

    # Add the atoms of the two solutes to the membrane structure
    # the coordinates of the atoms will be subsequently changed
    for i, atom in enumerate(solutefile1.atoms, 1) :
        atom.serial = len(boxfile.atoms)+i
        atom.residue = len(boxfile.residues)+1
    boxfile.extend(solutefile1)

    for i, atom in enumerate(solutefile2.atoms, 1) :
        atom.serial = len(boxfile.atoms)+i
        atom.residue = len(boxfile.residues)+1
    boxfile.extend(solutefile2)

    # Calculate the center of membrane
    memcoord = []
    for atom in boxfile.atoms :
        if atom.name.strip() == args.patom :
            memcoord.append(atom.xyz)
    memcent = np.asarray(memcoord).mean(axis=0)
    print "Liposome center = %.3f %.3f %.3f"%tuple(memcent)

    # Try to load a history of displacement vectors
    disp_vec = []
    have_old = True
    try :
        f = open("disp_vec.dat","r")
    except :
        have_old = False
    if have_old :
        disp_vec = [line.strip().split() for line in f.readlines()]
        f.close()
        disp_vec = np.array(disp_vec, dtype=float)

    # Find vector for insertion
    if have_old and args.nvec is not None :
        currvec = disp_vec[args.nvec-1,:]
    else :
        currvec = geo.sphere_surf_rand()

    # Write out the histoy of x,y positions
    with open("disp_vec.dat", "w") as f:
        if have_old :
            for vec in disp_vec :
                f.write("%8.3f %8.3f %8.3f\n"%tuple(vec))
        if not (have_old and args.nvec is not None) :
            f.write("%8.3f %8.3f %8.3f\n"%tuple(currvec))



    for ri, (rad1, rad2) in enumerate(rads_pack) :
        disp1 = memcent - solcent + rad1*currvec
        disp2 = memcent - solcent - rad2*currvec

        new1 = memcent + rad1*currvec
        new2 = memcent - rad2*currvec
        print "Put solute centers at: %.3f %.3f %.3f and %.3f %.3f %.3f"% \
            (new1[0], new1[1], new1[2], new2[0], new2[1], new2[2])

        for xyz, atom1, atom2 in zip(solxyz, solutefile1.atoms, solutefile2.atoms) :
            atom1.set_xyz(xyz + disp1)
            atom2.set_xyz(xyz + disp2)

        boxfile.write_gro("%s_r%d.gro"%(args.out,ri))
