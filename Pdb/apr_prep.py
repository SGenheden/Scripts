# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to prep a PDB structure for APR simulations
1.) Align the long axis of the host with the z-axis
2.) Align the short axis of the host with the y-axis
3.) Set origin of the G1 atom
4.) Add dummy atoms to the structure

The -a flag specifies the axes of the host
It takes to arguments and each argument is a list of atom names of the host.
The long axis goes from the center of the atoms in the two lists,
from the first centroid to the second centroid.
The short axis goes between the first two atoms in the first list

Examples
--------
  apr_prep.py prot.pdb -a C3,C29 C84,C78
  apr_prep.py prot.pdb -i apr.in -a C3,C29 C84,C78
"""

import argparse
import os
import sys

import numpy as np

from sgenlib import fitting
from sgenlib import geo
from sgenlib import pdb


def rotate_struct(struct, residue, axis, normal) :
    """
    Rotate the structure so that the axis of residue is aligned with normal
    """
    normv = np.zeros(3)
    normv[["x","y","z"].index(normal)] = 1.0

    # Move structure to center of coordinates
    res_xyz = residue.collect("xyz")
    center = res_xyz.mean(axis=0)
    res_xyz -= center
    xyz = struct.xyz - center
    xyz = xyz - center

    # Find the rotation matrix and rotate the full structure
    rotvec = geo.rotaxis(axis,normv)
    alpha = geo.angle(axis,normv)
    rotmat = geo.rotation_matrix(alpha,rotvec)
    xyz = fitting.rotate(xyz,rotmat)+center
    struct.update_xyz(xyz)

if __name__ == "__main__":


    # Setup a parser of the command-line arguments
    parser = argparse.ArgumentParser(description="Program to make prep a PDB file for APR simulations")
    parser.add_argument('file',help="the PDB file")
    parser.add_argument('-o','--out',help="the output filename")
    parser.add_argument('-i','--input',help="the APR input filename",default="apr.in")
    parser.add_argument('-a','--axis',nargs=2,help="the host axis definition")
    args = parser.parse_args()

    pdbfile = pdb.PDBFile(filename=args.file)

    with open(args.input, "r") as f :
        for line in f.readlines() :
            if line.strip().startswith("H1") :
                h1 = line.strip().split("=")[1].strip()
            elif line.strip().startswith("G1") :
                g1 = line.strip().split("=")[1].strip()
    g1 = pdbfile.atomindex(g1)
    hostname = h1.split('@')[0][1:]

    # Find host residue
    host = None
    for residue in pdbfile.residues :
        if residue.resname.strip() == hostname :
            host = residue
            break

    # Define the host axis
    host_xyz = host.collect("xyz")
    axis_indices = []
    centroids = []
    for ax in args.axis :
        coords = []
        indices = []
        for atm in ax.split(",") :
            idx = host.index_by_name(atm)
            indices.append(idx)
            coords.append(host_xyz[idx,:])
        axis_indices.append(indices)
        centroids.append(np.asarray(coords).mean(axis=0))

    # First rotate the structure around the long axis of the host
    rotate_struct(pdbfile, host, centroids[1]-centroids[0], 'z')

    # Then rotate around the short axis of the host
    host_xyz = host.collect("xyz")
    rotate_struct(pdbfile, host, host_xyz[axis_indices[0][0],:]-
            host_xyz[axis_indices[0][1],:], 'y')

    # Next move G1 to the new origin
    pdbfile.update_xyz(pdbfile.xyz-pdbfile.xyz[g1,:])

    # And finally add three new residues with dummy atoms
    new_xyz = [[0.000,3.500,-14.500],
                [0.000,0.000,-11.000],
                [0.000,0.000,-6.000]]
    for i in range(3) :
        res = pdb.Residue()
        atom = pdb.Atom()
        atom.idx = 1
        atom.hetatm = False
        atom.serial = 1
        atom.name = " Pb "
        atom.residue = 1
        atom.resname = "DUM"
        atom.set_xyz(new_xyz[i])
        res.append(atom)
        pdbfile.atoms.insert(0, atom)
        pdbfile.residues.insert(0, res)
    pdbfile.renumber()

    if args.out is None :
        args.out = os.path.splitext(args.file)[0]+"_apr.pdb"
    pdbfile.write(args.out)
