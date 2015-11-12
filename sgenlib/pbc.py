# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Routines to modify coordinates in periodic boxes.
Only works with rectangular geometries.
"""

import numpy as np

def unwrap_vector(dr,box):
    """
    Unwraps a vector
    """
    dr = np.divide(dr,box)
    dr = np.around(dr)
    return np.multiply(box,dr)

def make_whole(atoms,box) :
    """
    Make a set of Atom objects whole, i.e.
    make sure that they are not broken ower the periodic box

    Parameters
    ----------
    atoms : list of Atom objects
        the atoms to be modified
    box : numpy.ndarray
        the box dimensions
    """

    # Uwrap each atom, based on the position of the previous atom
    for atomi,atomj in zip(atoms[1:],atoms[:-1]) :
        dr = unwrap_vector(atomj.xyz-atomi.xyz,box)
        atomi.set_xyz(atomi.xyz+dr)

def make_whole_xyz(xyz,box,dim=[True,True,True]) :
    """
    Make a list of coordinates whole

    Parameters
    ----------
    xyz : numpy.ndarray
        the coordinates to be modified
    box : nump.ndarray
        the box dimensions
    dim : list of bool
        indicators to do unwrapping in only specific dimensions
    """

    dim = np.asarray(dim)
    xyz2 = np.array(xyz,copy=True)
    # Uwrap each coordinate, based on the position of the previous coordinate set
    for i in range(1,xyz.shape[0]) :
        dr = unwrap_vector(xyz2[i-1,:] - xyz[i,:], box)
        xyz2[i,dim] = xyz[i,dim] + dr[dim]
    return xyz2

def center(group1,group2,box) :
    """
    Center a set of coordinates by moving the secondary
    group of coordinates

    Parameters
    ----------
    group1 : list of Atom objects
        the atoms that should be central
    group2 : list of Residue objects
        the residues that should be moved
    box : nump.ndarray
        the box dimensions
    """
    com1 = np.asarray([atom.xyz for atom in group1]).mean(axis=0)
    for residue in group2 :
        com2 = np.asarray([atom.xyz for atom in residue.atoms]).mean(axis=0)
        dr = unwrap_vector(com1 - com2, box)
        for atom in residue.atoms :
            atom.set_xyz(atom.xyz+dr)

    # Make sure that the center of the coordinates is at the middle of the box
    delta = com1 - box/2.0
    for atom in group1 :
        atom.set_xyz(atom.xyz-delta)
    for residue in group2 :
        for atom in residue.atoms :
            atom.set_xyz(atom.xyz-delta)

def wrap(atoms,box) :
    """
    Wrap atom coordinates inside the box

    Experimental!

    Parameters
    ----------
    atoms : list of Atom objects
        the atoms to wrap
    box : nump.ndarray
        the box dimensions
    """

    DELTA = 0.001

    boxlow = np.zeros(3)
    boxhigh = np.zeros(3)

    box = np.asarray(box)
    if box.shape[0] == 3 :
        boxhigh = np.copy(box)
    elif box.shape[0] == 6 :
        boxlow = box[:3]
        boxhigh = box[3:]
    boxlen = boxhigh - boxlow

    for atom in atoms :
        xyz = np.copy(atom.xyz)
        for i in range(3) :
            if xyz[i] < boxlow[i]-DELTA :
                xyz[i] += boxlen[i]
            elif xyz[i] > boxhigh[i]+DELTA :
                xyz[i] -= boxlen[i]
        atom.set_xyz(xyz)
