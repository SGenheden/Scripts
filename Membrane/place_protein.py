# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to place a protein next to a membrane for "push"-simulation

Example:
  place_protein.py -prot prot.gro -mem membrane.gro
"""

import argparse
import os
import sys

import numpy as np

from sgenlib import pdb

def _clashes(xyz1, xyz2):

    for i in range(xyz1.shape[0]) :
        d = np.sum(np.power(xyz2-xyz1[i,:],2),axis=1)
        if d.min() < 2.0 :
            return True
    return False

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Program to place a protein next to a membrane")
    parser.add_argument('-prot','--protein',help="the protein PDB file")
    parser.add_argument('-mem','--membrane',help="the membrane PDB file")
    parser.add_argument('--lipidref',help="the name of the lipid reference atom",default="P")
    parser.add_argument('-o','--out',help="the output filename",default="placed.gro")
    args = parser.parse_args()

    prot = pdb.PDBFile(args.protein)
    mem = pdb.PDBFile(args.membrane)

    # Calculate the delta that will make the centroid of the two structures to overlap
    avmem = np.average([atom.xyz for atom in mem.atoms if atom.name.strip() == args.lipidref], axis=0)
    avprot = np.average(prot.xyz,axis=0)
    delta = avmem-avprot

    # Calculate how much from the edge of the membrane to shift the protein
    minprot = np.min(prot.xyz,axis=0)
    minatmx = np.argmin(prot.xyz, axis=0)[0]
    prot_lenx = (avprot-minprot)[0]
    mem_lenx = mem.xyz[:,0].max() - avmem[0]
    shift = np.array([prot_lenx+mem_lenx, 0, 0])

    # Picks out a few atoms for clashes check
    memsel = mem.xyz[:,0] > avmem[0] + 0.5 * mem_lenx
    protsel = prot.xyz[:,0] < avprot[0] - 0.5 * prot_lenx

    # Do the first shift
    prot.update_xyz(prot.xyz+delta+shift)
    # then shift it more in x until the protein clashes with the membrane
    shift = np.array([-2.0, 0.0, 0.0])
    while not _clashes(prot.xyz[protsel], mem.xyz[memsel]) :
        prot.update_xyz(prot.xyz+shift)

    # Calculate the new box size and centers the system in it
    extent_x = prot.xyz[:,0].max() - mem.xyz[:,0].min()
    newbox = np.array([extent_x, extent_x, mem.box[2]])
    newcenter = 0.5 * newbox
    avprot = np.average(prot.xyz,axis=0)
    avmem = np.average(mem.xyz,axis=0)
    currcenter = 0.5 * (avmem + avprot)
    delta = newcenter - currcenter
    delta[2] = 0.0 # Don't re-center in z-dimension, unnecessary
    mem.update_xyz(mem.xyz + delta)
    prot.update_xyz(prot.xyz + delta)

    prot.box = newbox
    prot.extend_residues(mem.residues)
    prot.renumber(doatoms=True, doresidues=True)
    prot.write(args.out)
