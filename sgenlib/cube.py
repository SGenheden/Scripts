# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Helper routines for dealing with CUBE files
"""

import numpy as np
import sys,os
import matplotlib.pylab as plt
import matplotlib.colors as colors
import matplotlib

#
# Read a file in CUBE format and return the grid, the atoms and the center of coordinates
#
def read_cube(filename) :

  f = open(filename,'r')
  # Skip two line
  f.readline()
  f.readline()
  # Third line is atom and center of box
  (natom,xcent,ycent,zcent) = f.readline().strip().split()
  natom = int(natom)
  cent = np.array([xcent,ycent,zcent],dtype=float)
  # Number of X bins and vector len
  (nx,xdim,dummy1,dummy2) = f.readline().strip().split()
  nx = int(nx)
  # Number of Y bins and vector len
  (ny,ydim,dummy1,dummy2) = f.readline().strip().split()
  ny = int(ny)
  # Number of Z bins and vector len
  (nz,zdim,dummy1,dummy2) = f.readline().strip().split()
  nz = int(nz)
  # Read atoms
  atoms = [f.readline().strip().split() for i in range(natom)]
  atoms = np.array(atoms,dtype=float)
  # Read the grid
  grid = np.zeros([nx,ny,nz])
  for x in range(nx) :
    for y in range(ny) :
      line = f.readline()
      while len(line) < 4 :
        line = f.readline()
      grid[x,y,:] = np.array(line.strip().split(),dtype=float)
  f.close()
  return (grid,atoms,cent)

#
# Read a file in CUBE format and manipulate it
#   1) crop the x and y dimensions to a maximum of 70 Angstromgs
#   2) sum the grid in z-dimensions
#   3) divide the sum into two halvs (upper and lower leaflet)#
#   4) convert atom coordinates to Angstroms
#   5) optional convert the density to a free energy (in kcal/mol)
#
def cube2mat(cubefile,matfile="",dgref=None) :

  grid,atoms,cent = read_cube(cubefile)
  atoms[:,2:] = atoms[:,2:]*0.529177249 # Convert to Angstroms

  MAT_MAX = 70
  xlow =  int(grid.shape[0]/2) - MAT_MAX / 2
  xhigh = int(grid.shape[0]/2) + MAT_MAX / 2
  ylow =  int(grid.shape[1]/2) - MAT_MAX / 2
  yhigh = int(grid.shape[1]/2) + MAT_MAX / 2
  midz  = int(grid.shape[2]/2) + 1

  grid = grid[xlow:xhigh,ylow:yhigh,:]

  low  = grid[:,:,:midz].sum(axis=2)
  upp  = grid[:,:,midz:].sum(axis=2)
  both = grid.sum(axis=2)

  if dgref != None :
    if dgref > 0 :
      low[low>0] = 1.987*0.3*np.log(low[low>0]/dgref)
      upp[upp>0] = 1.987*0.3*np.log(upp[upp>0]/dgref)
      both[both>0] = 1.987*0.3*np.log(both[both>0]/dgref)
    else :
      dgref = abs(dgref)
      low[low>0] = low[low>0]/dgref
      upp[upp>0] = upp[upp>0]/dgref
      both[both>0] = both[both>0]/dgref
  if len(matfile) > 0 :
    np.savez(matfile,atoms=atoms,cent=cent,lower=low,upper=upp,both=both)
  return {'atoms' : atoms,'cent' : cent,'lower' : low,'upper' : upp, 'both' : both}
