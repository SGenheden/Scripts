# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Routines to modify coordinates in periodic boxes.

Only works with rectangular geometries!
"""

import numpy as np
import os
import sys

thispath = os.path.dirname(os.path.abspath(__file__))
oneup = os.path.split(thispath)[0]
sys.path.insert(0,os.path.join(oneup,"Pdb"))
sys.path.insert(0,os.path.join(oneup,"Lammps"))

import pdb
import lammps

def make_whole(atoms,box,atom0idx=0,verbose=False) :
  """
  Make a set of Pdb or Lammps atom whole, i.e.
  make sure that they are not broken ower the periodic box

  Parameters
  ----------
  atoms : list of Atom objects 
    the atoms to be modified
  box : numpy.ndarray
    the box dimensions
  atom0idx : int, optional
    based the wrapping on this atom
    not used any more
  verbose : bool, optional
    turn on debug info
    not currently used
  """

  """
  Deprecated code

  ns = np.zeros([len(atoms),3])
  # Check for the atom which causes the least atoms to move
  if atom0idx < 0 :
    for j,atom0 in enumerate(atoms) : # Loop over all possible atom0
      for i,atomi in enumerate(atoms) :
        if i==j : continue
        dr = atom0.xyz-atomi.xyz
        dr = np.divide(dr,box)
        ns[j,:] = ns[j,:] + np.round(dr)
    
    # Check for the minimum move
    atom0idx = 0
    if not np.all(ns==0.0) :
      ns = ns.sum(axis=1)      
      atom0idx = np.argmin(ns)

  # Check periodicity relative to atom0
  atom0 = atoms[atom0idx]
  for i,atomi in enumerate(atoms) :
    if i==atom0idx : continue
    dr = atom0.xyz-atomi.xyz
    dr = np.divide(dr,box)
    dr = np.around(dr)
    dr = np.multiply(box,dr)
    atomi.x += dr[0]
    atomi.y += dr[1]
    atomi.z += dr[2]
    atomi.xyz = np.array([atomi.x,atomi.y,atomi.z])
    if verbose : print atomi"""

  # Uwrap each atom, based on the position of the previous atom
  for atomi,atomj in zip(atoms[1:],atoms[:-1]) :
    dr = atomj.xyz-atomi.xyz
    dr = np.divide(dr,box)
    dr = np.around(dr)
    dr = np.multiply(box,dr)
    atomi.x += dr[0]
    atomi.y += dr[1]
    atomi.z += dr[2]
    atomi.xyz = np.array([atomi.x,atomi.y,atomi.z])    

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
    dr = xyz2[i-1,:]-xyz[i,:]
    dr = np.divide(dr,box)
    dr = np.around(dr)
    dr = np.multiply(box,dr)
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
  group2 : Residue object
    the residues that should be moved
  box : nump.ndarray
    the box dimensions  
  """

  # Calculate the centroid of group1
  com1 = np.zeros(3)
  for atom in group1 :
    com1 = com1 + atom.xyz
  com1 = com1 / float(len(group1))
  
  # Move each residue as a unit
  for residue in group2 :
    com2 = np.zeros(3)
    for atom in residue.atoms :
      com2 = com2 + atom.xyz
    com2 = com2 / float(len(residue.atoms))
    
    dr = com1 - com2
    dr = np.divide(dr,box)
    dr = np.round(dr)
    dr = np.multiply(box,dr)
    for atom in residue.atoms :
      atom.x += dr[0]
      atom.y += dr[1]
      atom.z += dr[2]
      atom.xyz = np.array([atom.x,atom.y,atom.z])

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
        
