# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to insert a solute at different depths in a membrane

Two solutes will be inserted separated by half a box length in
the membrane plane and the total separation in z.

The solute will all likely overlap with the membrane and therefore
it is essential to grow it before attemping e.g. umbrella sampling

Example:
  insert_in_membrane.py -f aceticacid.pdb -b membrane.gro -z {0..37}
"""

import argparse

import numpy as np
import numpy.random as random

from sgenlib import pdb

if __name__ == '__main__' :

    parser = argparse.ArgumentParser(description="Insert a solute in membrane for umbrella sampling")
    parser.add_argument('-f','--file',help="the solute file")
    parser.add_argument('-b','--box',help="the membrane box")
    parser.add_argument('-o','--out',help="the output file",default="inserted")
    parser.add_argument('-z','--zvals',nargs="+",type=float,help="the position in the z-direction")
    parser.add_argument('-p','--patom',help="the name of the phosphate atom")
    parser.add_argument('-n','--nxy',type=int,help="the index of a previously employed xy-position")
    args = parser.parse_args()

    zvals_rev = args.zvals[::-1]
    zvals_pack = zip(args.zvals,zvals_rev)

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
    if args.patom is None :
        args.patom = "P"
        if len([atom for atom in boxfile.atoms
                if atom.name.strip() == args.patom]) == 0 :
            args.patom = "PO4"
    memz = []
    for atom in boxfile.atoms :
        if atom.name.strip() == args.patom :
            memz.append(atom.z)
    memcent = np.asarray(memz).mean()
    print "Membrane center = %.3f"%memcent

    # Try to load a history of x,y positions
    xy_pos = []
    have_old = True
    try :
        f = open("xy_pos.dat","r")
    except :
        have_old = False
    if have_old :
        xy_pos = [line.strip().split() for line in f.readlines()]
        f.close()
        xy_pos = np.array(xy_pos, dtype=float)

    # Find xy position to insert
    if have_old and args.nxy is not None :
        x1 = xy_pos[args.nxy-1, 0]
        y1 = xy_pos[args.nxy-1, 1]
    else :
        x1 = (boxfile.box[0]-10.0) * random.random_sample() + 5.0
        y1 = (boxfile.box[1]-10.0) * random.random_sample() + 5.0
        if have_old :
            diff = (xy_pos - np.array([x1,y1]))**2
            while np.any(diff<25.0) :
                x1 = (boxfile.box[0]-10.0) * random.random_sample() + 5.0
                y1 = (boxfile.box[1]-10.0) * random.random_sample() + 5.0
                diff = (xy_pos - np.array([x1,y1]))**2

    # Write out the histoy of x,y positions
    with open("xy_pos.dat", "w") as f:
        if have_old :
            for pos in xy_pos :
                f.write("%8.3f %8.3f\n"%(pos[0],pos[1]))
        if not (have_old and args.nxy is not None) :
            f.write("%8.3f %8.3f\n"%(x1,y1))


    if x1 > boxfile.box[0] / 2.0 :
        x2 = x1 - boxfile.box[0] / 2.0
    else :
        x2 = x1 + boxfile.box[0] / 2.0

    if y1 > boxfile.box[1] / 2.0 :
        y2 = y1 - boxfile.box[1] / 2.0
    else :
        y2 = y1 + boxfile.box[1] / 2.0

    for zval1, zval2 in zvals_pack :
        z1 = memcent + zval1
        disp1 = np.array([x1, y1, z1]) - solcent

        z2 = memcent - zval2
        disp2 = np.array([x2, y2, z2]) - solcent

        print "Put solute centers at: %.3f %.3f %.3f and %.3f %.3f %.3f"% \
            (x1, y1, z1, x2, y2, z2)

        for xyz, atom1, atom2 in zip(solxyz, solutefile1.atoms, solutefile2.atoms) :
            atom1.set_xyz(xyz + disp1)
            atom2.set_xyz(xyz + disp2)

        boxfile.write_gro("%s_z%d.gro"%(args.out,zval1))
