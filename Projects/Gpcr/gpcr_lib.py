# Author: Samuel Genheden samuel.genheden@gmail.com

"""
This module contain classes and routines
to analyse GPCR simulations
"""

import os
import sys
import struct
from ConfigParser import SafeConfigParser

import numpy as np
from scipy.spatial.distance import cdist
from scipy.ndimage.measurements import center_of_mass
import MDAnalysis.lib.distances as distances
from Bio.KDTree import KDTree
import matplotlib as plt

# Import the PDB and Plot modules
from sgenlib import pdb
from sgenlib import colors

# Change this to load protein template and Xray structures from a different location
PROT_INFO_PATH = "/Users/samuel/Dropbox/Research/GPCR/Prot_info/"

# Labels for the lower and upper densities
side_name = {"low":"intra.","upp":"extra."}

class GridTemplate(object) :
  """
  Class to store the template of 2D grid

  Attributes
  ----------
  resolution : float
    the resolution, i.e. the size of each pixel
  size : tuple of integers
    the size of the grid
  edgesx : NumpyArray
    the edges along the x-dimension
  edgesy : NumpyArray
    the edges along the y-dimension
  origin : tuple of float
    the origin of the grid
  len : tuple of float
    the length in real dimensions
  """
  def __init__(self,xyz,resolution=1.0) :
    coord_min = xyz[:,:2].min(axis=0)
    coord_max = xyz[:,:2].max(axis=0)
    coord_len = coord_max - coord_min
    self.size = tuple(np.ceil(coord_len/resolution).astype(np.int_))
    self.edgesx = np.linspace(coord_min[0],coord_max[0],self.size[0],endpoint=True)
    self.edgesy = np.linspace(coord_min[1],coord_max[1],self.size[1],endpoint=True)
    self.size = tuple(np.ceil(coord_len/resolution).astype(np.int_)+1)
    self.origin = tuple(coord_min)
    self.len = tuple(coord_len)
    self.resolution = resolution
  def __str__(self) :
    str = "Resolution = %.1f A\n"%self.resolution
    str += "Size = %d x %d\n"%self.size
    str += "Edges x-limits = %.3f %.3f\n"%(self.edgesx[0],self.edgesx[-1])
    str += "Edges y-limits = %.3f %.3f\n"%(self.edgesy[0],self.edgesy[-1])
    str += "Origin = %.3f , %.3f\n"%self.origin
    str += "Length = %.3f , %.3f\n"%self.len
    return str
  def create(self,depth=None) :
    """
    Create a numpy array of appropriate size
    """
    if depth is None :
      return np.zeros(self.size)
    else :
      return np.zeros([self.size[0],self.size[1],depth])
  def indices(self,xyz) :
    """
    Returns the grid coordinates for a set of Cartesian coordinates
    """
    xidx = np.digitize(xyz[:,0],self.edgesx)
    yidx = np.digitize(xyz[:,1],self.edgesy)
    return np.array([xidx,yidx])
  def save(self,mats,size,filename) :
    """
    Save a dictionary of matrices to a npz-file

    Truncates them to a specific size
    """
    mats2 = {}
    for m in mats :
      mats2[m] = np.array(mats[m],copy=True)
      if mats2[m].shape == self.size  :
        mats2[m] = self.truncate(mats2[m],size)
    np.savez(filename,**mats2)
  def truncate(self,mat,size) :
    """
    Truncate a numpy array to a specific size around the centre
    """
    centx = int(float(self.size[0])/2.0)
    centy = int(float(self.size[1])/2.0)
    lowx = centx-size
    lowy = centy-size
    highx = centx+size
    highy = centy+size
    return mat[lowx:highx,lowy:highy]

class ProteinTemplate(object) :
  """
  Class to store residue and atom
  template information for a protein.
  This include translation from internal indices to Xray indices

  Attributes
  ----------
  residues : list of integers
    the residue number in the reference
  mob2ref : dictionary
    look-up from mobile to reference number
  codes : list of strings
    the residue codes
  names : list of strings
    the residue names
  helices : list
    the first and last atom indices for residue ranges (e.g. 7 transmembrane helices)
  rhelices : list
    the first and last residue indices for residue ranges
  """
  def __init__(self,filename) :

    parser = SafeConfigParser()
    parser.read(filename)

    self.residues = map(int,parser.get("Residue match","Numbers-ref").split())
    self.mob2ref = {}
    for r in self.residues :
      v = int(parser.get("Residue match","%d"%r))
      self.mob2ref[v] = r

    self.codes = parser.get("Residue names","Codes").split()
    self.names = parser.get("Residue names","Names").split()

    deflist = parser.get("Residue ranges","Def-list")[1:-1]
    deflist = deflist.replace(", ",":")
    self.rhelices = []
    for d in deflist.split(",") :
      r1,r2 = d[1:-1].split(":")
      self.rhelices.append([int(r1),int(r2)])

    defalist = parser.get("Residue ranges","Def-atom-list")[1:-1]
    defalist = defalist.replace(", ",":")
    self.helices = []
    for d in defalist.split(",") :
      a1,a2 = d[1:-1].split(":")
      self.helices.append([int(a1),int(a2)])

####################
# Analysis routines
####################

def anal_contacts(prot,rfirst,mols,cutoff,midz,box,molfile,resfile,
                    burresfile=None,buried=None,reslist=None,reslistfile=None) :
  """
  Calculates a range of contacts and write out contact vectors to files

  Parameters
  ----------
  prot : AtomSelection
    the protein atoms
  rfirst : NumpyArray
    the first atom of each residue in the protein
  mols : NumpyArray
    the centroid of the molecule of interest
  cutoff : float
    the contact cut-off
  midz : float
    the middle of the bilayer
  box : NumpyArray
    the box sides
  molfile : fileobject
    the file to write molecular contacts to
  resfile : fileobject
    the file to write residue contacts to
  burresfile : fileobject, optional
    the file to write buried residue contacts to
  buried : NumpyArray of boolean, optional
    flag to indicate buried molecule
  reslist : list
    a list of residues to write out individual mol contacts
  reslistfile : fileobject
    the file to write out individual mol contacts to
  """
  molon = np.zeros(len(mols),dtype=bool)
  reson = np.zeros(len(rfirst)-1,dtype=bool)
  if buried is not None :
    burreson = np.zeros(len(rfirst)-1,dtype=bool)
  if reslist is not None :
    resonlist = np.zeros([len(reslist),len(mols)],dtype=bool)

  # Calculate all distances at once
  if box is not None :
    dist_all = distances.distance_array(np.array(mols), prot.get_positions(), box)
  else :
    kdtree = KDTree(dim=3, bucket_size=10)
    kdtree.set_coords(prot.get_positions())

  # Loop over all molecules
  for mi,mol in enumerate(mols) :
    if box is not None :
      dist = dist_all[mi,:]
    else :
      dist = np.ones(prot.get_positions().shape[0])+cutoff
      kdtree.search(np.array([mol]), cutoff)
      for i in kdtree.get_indices() : dist[i] = 0.0

    # Check if this molecule is on
    molon[mi] = dist.min() < cutoff
    if  molon[mi] :
      # Check contacts for reach residue
      for i in range(len(rfirst)-1) :
        if dist[rfirst[i]:rfirst[i+1]].min() < cutoff :
          if (burresfile is not None and buried[mi]) :
            burreson[i] = True
          else :
            reson[i] = True
          if reslist is not None and i in reslist:
            resonlist[reslist.index(i),mi] = True

  # Write state information to file
  write_booleans(molfile,molon)
  write_booleans(resfile,reson)
  if buried is not None :
    write_booleans(burresfile,burreson)
  if reslist is not None:
    write_booleans(reslistfile,resonlist.reshape(len(reslist)*len(mols)))

def anal_density(low_sel,count,grid,mat,buried=None,**kwargs) :
  """
  Put the coordinates of molecules on a grid

  By default, two grids are accumulated; one for the molecules in the lower leaflet
  and one for the molecules in the upper leaflet
  If buried is given, the split of the molecules is made three-fold, in addition to
  molecules in the lower and upper leaflet, molecules in the middle of the bilayer
  is put on a separate grid. These molecules are not put on any other grid.

  Parameters
  ----------
  low_sel : NumpyArray of booleans
    selection of molecules in the lower leaflet
  count : string
    the matrix element that holds the total count of discretised molecules
  grid : GridTemplate
    the grid definition, used for discretisation
  mat : dictionary of NumpyArray
    the accumulated grids for each molecule
    this matrix will be modified by this routine
  buried : NumpyArray of booleans, optional
    selection of molecules in the bilayer middle
  **kwargs : dictionary of NumpyArrays
    can be either a discretised coordinates or Cartesian coordinates
  """

  # Select molecule in the upper leaflet as the negative of the selection in
  # lower leaflet
  upp_sel = np.logical_not(low_sel)
  if buried is not None :
    # Remove buried molecules from lower and upper selections
    low_sel = np.logical_and(low_sel,np.logical_not(buried))
    upp_sel = np.logical_and(upp_sel,np.logical_not(buried))
  # Accumulate total number of discretised molecules
    mat["cnt_"+count][2] = mat["cnt_"+count][2] + buried.sum()
  mat["cnt_"+count][0] = mat["cnt_"+count][0] + low_sel.sum()
  mat["cnt_"+count][1] = mat["cnt_"+count][1] + upp_sel.sum()

  # Now accumulate grids for each molecule
  for k in kwargs :
    if "low_"+k not in mat : continue
    # v can either be a 2 x n array in case it is assumed to be discretised coordinates
    # if v is a n x 3 array it is assumed to be Cartesian coordinates
    v = kwargs[k]
    # If it is Cartesian coordinates, we need to discretise it first
    if v.shape[1] == 3 :
      v = grid.indices(v)
    lidx = v[:,low_sel]
    uidx = v[:,upp_sel]
    mat["low_"+k][lidx[0],lidx[1]] = mat["low_"+k][lidx[0],lidx[1]] + 1
    mat["upp_"+k][uidx[0],uidx[1]] = mat["upp_"+k][uidx[0],uidx[1]] + 1
    if buried is not None :
      midx = v[:,buried]
      mat["mid_"+k][midx[0],midx[1]] = mat["mid_"+k][midx[0],midx[1]] + 1

def anal_jointdist(prot,rfirst,mols,cutoff,midz,box,mat,matburr,buried=None) :
  """
  Calculates the joint probability of residue interactions

  Parameters
  ----------
  prot : AtomSelection
    the protein atoms
  rfirst : NumpyArray
    the first atom of each residue in the protein
  mols : NumpyArray
    the centroid of the molecule of interest
  cutoff : float
    the contact cut-off
  midz : float
    the middle of the bilayer
  box : NumpyArray
    the box sides
  mat : NumpyArray
    the joint probability for non-buried molecules
  matburr : NumpyArray
    the joint probability for buried molecules
  buried : NumpyArray of boolean, optional
    flag to indicate buried molecule
  """
  imat = np.zeros(mat.shape)
  imatburr = np.zeros(matburr.shape)

  # Calculate all distances at once
  if box is not None :
    dist_all = distances.distance_array(np.array(mols), prot.get_positions(), box)
  else :
    kdtree = KDTree(dim=3, bucket_size=10)
    kdtree.set_coords(prot.get_positions())

  # Loop over all molecules
  for mi,mol in enumerate(mols) :
    if box is not None :
      dist = dist_all[mi,:]
    else :
      dist = np.ones(prot.get_positions().shape[0])+cutoff
      kdtree.search(np.array([mol]), cutoff)
      for i in kdtree.get_indices() : dist[i] = 0.0

    # Check if this molecule is on
    if dist.min() >= cutoff : continue

    # Check contacts for reach residue
    for i in range(len(rfirst)-1) :
      if dist[rfirst[i]:rfirst[i+1]].min() >= cutoff : continue
      if (buried is not None and buried[mi]) :
        imatburr[i,i] = 1
      else :
        imat[i,i] = 1

      for j in range(i+1,len(rfirst)-1) :
        if dist[rfirst[j]:rfirst[j+1]].min() >= cutoff : continue

        if (buried is not None and buried[mi]) :
          imatburr[i,j] = 1
          imatburr[j,i] = 1
        else :
          imat[i,j] = 1
          imat[j,i] = 1

  return (mat+imat,matburr+imatburr)

def anal_orderparams(heads,head_idx,midz,bonds,mat) :
  """
  Computes the tail order parameters

  Parameters
  ----------
  heads : AtomSelection
    the head particles
  head_idx : NumpyArray
    the discretised coordinates of the head beads
  midz : float
    the middle of the bilayer
  bonds : list of AtomSelection
    the bond selections
  mat : dictionary of NumpyArray
    the accumulate grids
    this will be modified
  """
  def calc_orderparam(atoms1,atoms2,norm) :
    """
    Calculate the order parameters for all lipids given a bond
    """
    vec = atoms1 - atoms2
    proj = np.multiply(vec,norm).sum(axis=1)**2 / np.sum(vec**2,axis=1)
    return proj

  bilayer_norm = np.array([0.0,0.0,1.0])
  low_sel = heads.get_positions()[:,2] < midz
  upp_sel = np.logical_not(low_sel)
  lidx = head_idx[:,low_sel]
  uidx = head_idx[:,upp_sel]

  for bi,bond in enumerate(bonds) :
    atoms1 = bond[0].get_positions()
    atoms2 = bond[1].get_positions()
    mat["all_sum"][bi,0] = mat["all_sum"][bi,0] + 0.5*(3.0*calc_orderparam(atoms1,atoms2,bilayer_norm).mean()-1)
    mat["all_sum"][bi,1] = mat["all_sum"][bi,1] + 0.5*(3.0*calc_orderparam(atoms1[low_sel,:],atoms2[low_sel,:],bilayer_norm).mean()-1)
    mat["all_sum"][bi,2] = mat["all_sum"][bi,2] + 0.5*(3.0*calc_orderparam(atoms1[upp_sel,:],atoms2[upp_sel,:],bilayer_norm).mean()-1)

    # Calculate the order parameters of the bond in the lower leaflet
    mat["low_sum"][lidx[0],lidx[1],bi]= mat["low_sum"][lidx[0],lidx[1],bi] + calc_orderparam(atoms1[low_sel,:],atoms2[low_sel,:],bilayer_norm)
    mat["low_count"][lidx[0],lidx[1],bi]  = mat["low_count"][lidx[0],lidx[1],bi]  + 1

    # Calculate the order parameter of the bond in the upper leaflet
    mat["upp_sum"][uidx[0],uidx[1],bi]= mat["upp_sum"][uidx[0],uidx[1],bi] + calc_orderparam(atoms1[upp_sel,:],atoms2[upp_sel,:],bilayer_norm)
    mat["upp_count"][uidx[0],uidx[1],bi] = mat["upp_count"][uidx[0],uidx[1],bi]  + 1

def anal_savemols(selection,nsnap,crd) :
  """
  Write out a selection of residues in xyz format

  Parameters
  ----------
  selection : AtomGroup
    the selected residues
  nsnap : integer
    the current snapshot
  crd : dictionary
    "frame" is the total number of residues printed out
    "file" is the fileobject to write to
  """
  for residue in selection.residues () :
    crd["frame"] += 1
    crd["file"].write("%d\n%d:%d\n"%(len(residue),mdcrd["frame"],nsnap))
    for atom in residue :
      crd["file"].write("%s %.3f %.3f %.3f\n"%(atom.name,atom.position[0],atom.position[1],atom.position[2]))

def anal_thickness(heads,head_idx,box,midz,mat) :
  """
  Computes the bilayer thickness and discretise it on a grid

  Parameters
  ----------
  heads : AtomSelection
    the head particles
  head_idx : NumpyArray
    the discretised coordinates of the head beads
  midz : float
    the middle of the bilayer
  mat : dictionary of NumpyArray
    the accumulate grids
    this will be modified
  """
  def calc_zdist(coords1,coords2) :
    """
    Calculate the z-distance between all lipids in one leaflet and the closest lipid in the other leaflet
    """
    dist = cdist(coords1[:,:2],coords2[:,:2],'sqeuclidean')
    j = np.argmin(dist,axis=1)
    return np.sqrt((coords2[j,2]-coords1[:,2])**2)

  low_sel = heads.get_positions()[:,2] < midz
  upp_sel = np.logical_not(low_sel)
  lidx = head_idx[:,low_sel]
  uidx = head_idx[:,upp_sel]
  low_coords = heads.get_positions()[low_sel,:]
  upp_coords = heads.get_positions()[upp_sel,:]

  # Calculate the thickness from the point of view of the lower leaflet
  zdist = calc_zdist(low_coords,upp_coords)
  mat["lower"][lidx[0],lidx[1]] = np.minimum(zdist,mat["lower"][lidx[0],lidx[1]])

  # Calculate the thickness from the point of view of the upper leaflet
  zdist = calc_zdist(upp_coords,low_coords)
  mat["upper"][uidx[0],uidx[1]] = np.minimum(zdist,mat["upper"][uidx[0],uidx[1]])

####################
# Density routines
####################

class Density :
  """
  Class to store a 2D molecular number density

  Attributes
  ----------
  mat : NxXxY Numpy array
    N individual densities with size X and Y
  nmat : integer
    the number of individual densities, set at construction
  excluded : XxY Numpy array
    voxels that are never occupied
  av : XxY Numpy array
    average density
  std : XxY Numpy array
    standard deviation of density
  dg : XxY Numpy array
    the free energy density
  dge : XxY Numpy array
    the standard deviation of the free energy
  nsnap : Numpy array
    number of snapshots for each density matrix
  molcount : Numpy array
    number of molecules forms the basis for each density
  """
  def __init__(self,mat,nfiles,**kwargs) :
    self.mat = np.zeros([nfiles,mat.shape[0],mat.shape[1]])
    self.mat[0,:,:] = mat
    self.nmat = 1
    self.excluded = np.ones(mat.shape)
    self.excluded[mat>0] = 0.0
    self._excluded_im = None
    self.av = None
    self.std = None
    self.dg = None
    self.dge = None
    if "nsnap" in kwargs :
      self.nsnap = np.zeros(nfiles)
      self.nsnap[0] = kwargs["nsnap"]
    else :
      self.nsnap = 1
    if "molcount" in kwargs :
      self.molcount = np.zeros(nfiles)
      self.molcount[0] = kwargs["molcount"]
    else :
      self.molcount = None
  def append(self,mat,**kwargs) :
    """
    Add another density to this collection

    Returns if the collection is already full
    """
    if self.nmat == self.mat.shape[0] : return
    self.mat[self.nmat,:,:] = mat
    if self.nsnap is not None and "nsnap" in kwargs :
      self.nsnap[self.nmat] = kwargs["nsnap"]
    if self.molcount is not None and "molcount" in kwargs :
      self.molcount[self.nmat] = kwargs["molcount"]
    self.nmat = self.nmat + 1
    self.excluded[mat>0] = 0.0
  def average(self) :
    """ Average the density collection
    """
    self.av = np.average(self.mat,axis=0)
    self.std = np.std(self.mat,axis=0)/np.sqrt(self.mat.shape[0])
  def calc_gibbs(self,cutoff=None,temp=310.0) :
    """
    Calculates the free energy density in kJ/mol

    Calls average() if it has been called previously

    Attributes
    ----------
    cutoff : float, optional
      the density cut-off, any density lower than this will be re-set to this value
    temp : float, optional
      the temperature in Kelvin
    """
    if self.av is None :
      self.average()
    factor = -1.987*temp/1000.0*4.184
    self.dg = np.zeros(self.av.shape)
    self.dge = np.zeros(self.av.shape)
    if cutoff is None :
      cutoff = 0.0
    else :
      sel = self.av < cutoff
      self.dg[sel] = factor*np.log(cutoff)
    sel = self.av > cutoff
    self.dg[sel] = factor*np.log(self.av[sel])
    self.dge[sel] = factor*np.divide(self.std[sel],self.av[sel])
  def centroid(self, rotate2D=None):
    if rotate2D is None:
        return center_of_mass(self.excluded)
    else:
        return center_of_mass(rotate2D(self.excluded))
  def cutoff_av(self,cutoff) :
    """
    Set average density below a cut-off to that cut-off value
    """
    sel = self.av < cutoff
    self.av[sel] = cutoff
  def erase_low(self,level) :
    """
    Set alpha level to 1.0 (opaque) for any voxel very density is low,
    as this will be rendered as white, it will "erase" low density

    Attributes
    ----------
    level : float
      the cut-off value of density
    """
    self.excluded_im()
    self._excluded_im[self.av<level,3] = 1.0
  def excluded_im(self) :
    """
    Creates white image for excluded density voxels, and sets alpha values to
    1 (opaque). However, for voxels that has density the alpha value is set to 0.0 (transparent)
    """
    if self._excluded_im is None :
      self._excluded_im = np.ones([self.excluded.shape[0],self.excluded.shape[1],4])
      self._excluded_im[self.excluded==0.0,3] = 0.0
    return self._excluded_im
  def _get_mat(self,mat) :
    """
    Returns the copy of a density matrix
    """
    if mat == "average" :
      return np.array(self.av,copy=True)
    elif mat == "std" :
      return np.array(self.std,copy=True)
    elif mat == "dg" :
      return np.array(self.dg,copy=True)
    elif mat == "dge" :
      return np.array(self.dge,copy=True)
  def max(self,mat) :
    """
    Returns the maximum of a specific density matrix, excluding zero density voxels
    """
    if self.av is None : self.average()
    density = self._get_mat(mat)
    return density[density!=0.0].max()
  def min(self,mat) :
    """
    Returns the minimum of a specific density matrix, excluding zero density voxels
    """
    if self.av is None : self.average()
    density = self._get_mat(mat)
    return density[density!=0.0].min()
  def plot(self,mat,axis,cmap,minval,maxval=None,reverseY=True, rotate2D=None) :
    """
    Plot the density

    Attributes
    ----------
    mat : string
      the density to plot
    axis : Axis object
      the axis to do the plotting on
    cmap : Cmap object
      the color map to use for the density
    minval : float
      the minimum value of the density to plot
    maxval : float, optional
      the maximum value of the density to plot
    reverseY : bool, optional
      whether to flip Y axis
    rotate2D : func, optional
      function to rotate the density in 2D
    Returns
    -------
    the image instance returned by imshow
    """
    density = self._get_mat(mat)
    if rotate2D is not None : density = rotate2D(density)
    if reverseY : density = density[:,::-1]
    if maxval is None :
      im =  axis.imshow(density,cmap=cmap,vmin=minval,extent=(-35,35,-35,35),origin="lower")
    else :
      im = axis.imshow(density,cmap=cmap,vmin=minval,vmax=maxval,extent=(-35,35,-35,35),origin="lower")
    density = np.array(self.excluded_im(),copy=True)
    if rotate2D is not None : density = rotate2D(density)
    if reverseY : density = density[:,::-1]
    axis.imshow(density,extent=(-35,35,-35,35),origin="lower")
    return im
  def scale(self,factor=None) :
    """
    Scale all densities with a factor
    """
    if factor is None :
      if self.mat is None : return
      for i in range(self.nmat) :
        self.mat[i,:,:] = self.mat[i,:,:] / self.nsnap[i]
    else :
      self.mat = self.mat / float(factor)
  def scale_by_uniform(self,nmol=None,boxlen=None,multifac=1) :
    """
    Scale all densities by uniform density in a box, excluding non-occupying voxels
    """
    if boxlen is None : boxlen = float(self.mat.shape[1])
    if nmol is None :
      if self.molcount is None : return
      uniform_density = self.molcount / float(boxlen*boxlen - self.excluded.sum())
      for i in range(self.nmat) :
        self.mat[i,:,:] = self.mat[i,:,:] / uniform_density[i]
    else :
      uniform_density = float(nmol)/float(boxlen*boxlen - self.excluded.sum())
      self.mat = self.mat / (multifac * uniform_density)
  def subtract(self,mat,tosubtract) :
    """ Subtract the density from another density
    """
    density1 = self._get_mat(mat)
    density2 = tosubtract._get_mat(mat)
    self.dg = density1 - density2
  def crop(self,xlen,ylen) :
    """ Crop the density to specific shape
    """
    if xlen == self.mat.shape[1] and ylen == self.mat.shape[2] : return
    centx = int(float(self.mat.shape[1])/2.0)
    centy = int(float(self.mat.shape[2])/2.0)
    lowx = centx-xlen/2
    highx = centx+xlen/2
    lowy = centy-ylen/2
    highy = centy+ylen/2
    self.mat = self.mat[:,lowx:highx,lowy:highy]
    self.excluded = self.excluded[lowx:highx,lowy:highy]
    if not self._excluded_im is None :
      self._excluded_im = self._excluded_im[lowx:highx,lowy:highy,:]
    if not self.av is None :
      self.av = self.av[lowx:highx,lowy:highy]
      self.std = self.std[lowx:highx,lowy:highy]
    if not self.dg is None :
      self.dg = self.dg[lowx:highx,lowy:highy]
      self.dge = self.dge[lowx:highx,lowy:highy]
  def export(self,mat) :
    return self._get_mat(mat),np.array(self.excluded_im(),copy=True)

class Densities :
  """
  Class to store a density for lower, upper and middle leaflet

  Attributes
  ----------
  low, upp, mid : Density
    the lower, upper and middle density
  lownam, uppnam, midnam : string
    the name of the density
  countnam : string
    the name of the count array
  """
  def __init__(self,filename,nfiles,matnam,countnam) :
    matdict = np.load(filename)

    self.lownam = "low_"+matnam
    self.uppnam = "upp_"+matnam
    self.midnam = "mid_"+matnam
    self.cntnam = "cnt_"+countnam

    nsnap = float(matdict["nsnap"])
    self.low = Density(matdict[self.lownam],nfiles,nsnap=nsnap,molcount=matdict[self.cntnam][0] / nsnap)
    self.upp = Density(matdict[self.uppnam],nfiles,nsnap=nsnap,molcount=matdict[self.cntnam][1] / nsnap)
    if matdict[self.cntnam].shape[0] > 2 :
      self.mid = Density(matdict[self.midnam],nfiles,nsnap=nsnap,molcount=matdict[self.cntnam][2] / nsnap)
    else :
      self.mid = None
      self.midnam = None
  def __getitem__(self,key) :
    if key.lower()[0:3] == "low" :
      return self.low
    elif key.lower()[0:3] == "upp" :
      return self.upp
    elif key.lower()[0:3] == "mid" :
      return self.mid
    else :
      raise Exception("Unknown density key: %s"%key)
  def average(self) :
    self.low.average()
    self.upp.average()
    if self.mid is not None : self.mid.average()
  def calc_gibbs(self,cutoff=None,temp=310.0) :
    self.low.calc_gibbs(cutoff=cutoff,temp=temp)
    self.upp.calc_gibbs(cutoff=cutoff,temp=temp)
    if self.mid is not None : self.mid.calc_gibbs(cutoff=cutoff,temp=temp)
  def cutoff_av(self,cutoff) :
    self.low.cutoff_av(cutoff)
    self.upp.cutoff_av(cutoff)
    if self.mid is not None : self.mid.cutoff_av(cutoff)
  def erase_low(self,level=1E-2) :
    self.low.erase_low(level)
    self.upp.erase_low(level)
    if self.mid is not None : self.mid.erase_low(level)
  def max(self,mat) :
    return max(self.low.max(mat),self.upp.max(mat))
  def min(self,mat) :
    return min(self.low.min(mat),self.upp.min(mat))
  def read(self,filename) :
    matdict = np.load(filename)
    nsnap = float(matdict["nsnap"])
    self.low.append(matdict[self.lownam],nsnap=nsnap,molcount=matdict[self.cntnam][0] / nsnap)
    self.upp.append(matdict[self.uppnam],nsnap=nsnap,molcount=matdict[self.cntnam][0] / nsnap)
    if self.mid is not None : self.mid.append(matdict[self.midnam],nsnap=nsnap,molcount=matdict[self.cntnam][2] / nsnap)
  def scale(self,factor=None) :
    self.low.scale(factor)
    self.upp.scale(factor)
    if self.mid is not None : self.mid.scale(factor)
  def scale_by_uniform(self,boxlen=None,nmol=None,ratio=0.5) :
    if nmol is None :
      self.low.scale_by_uniform(boxlen=boxlen)
      self.upp.scale_by_uniform(boxlen=boxlen)
      if self.mid is not None : self.mid.scale_by_uniform(boxlen=boxlen,multifac=2.0)
    else :
      self.low.scale_by_uniform(nmol=nmol*ratio,boxlen=boxlen)
      self.upp.scale_by_uniform(nmol=nmol*(1-ratio),boxlen=boxlen)
  def subtract(self,mat,tosubtract) :
    self.low.subtract(mat,tosubtract.low)
    self.upp.subtract(mat,tosubtract.upp)
    if self.mid is not None and tosubtract.mid is not None : self.mid.subtract(mat,tosubtract.mid)
  def crop(self,xlen,ylen) :
    self.low.crop(xlen,ylen)
    self.upp.crop(xlen,ylen)
    if self.mid is not None : self.mid.crop(xlen,ylen)
  def write(self,mat,filename) :
    lowmat,lowfilter = self.low.export(mat)
    uppmat,uppfilter = self.upp.export(mat)
    np.savez(filename,lowmat=lowmat,lowfilter=lowfilter,uppmat=uppmat,uppfilter=uppfilter)

def standard_read_and_process(filenames,densnam,countnam=None) :
  """
  Read a number of densities from file and process them in a standardized way
  """
  if countnam is None :
    countnam = densnam
  densities = Densities(filenames[0],len(filenames),densnam,countnam)
  for filename in filenames[1:] :
    densities.read(filename)
  densities.scale()
  densities.scale_by_uniform()
  densities.average()
  densities.erase_low()
  densities.calc_gibbs(cutoff=0.15)
  densities.crop(70,70)
  return densities

####################
# Plotting routines
####################

class XrayDensity :
  """
  Class to draw an X-ray structure on a grid
  """
  def __init__(self,filename,template,sigmafile=None) :

    # Read in a structure from file
    cholxyz = []
    protxyz = []
    protflag = []
    self.pdbfile = pdb.PDBFile(filename)
    flag = {"BB" : 1, "SC" :-1}
    for atom in self.pdbfile.atoms :
      if atom.resname[0:3] == "CHO" :
        cholxyz.append(atom.xyz)
      elif atom.name.strip()[:2] in ["BB","SC"] :
        protflag.append(flag[atom.name.strip()[:2]])
        protxyz.append(atom.xyz)
    self.box = self.pdbfile.box
    self.cholxyz = np.array(cholxyz)
    self.protxyz = np.array(protxyz)
    self.protflag  = np.array(protflag,dtype=int)
    self.helices  = template.helices
    self.template = template

    self.cholleaf = ""
    if self.cholxyz.shape[0] > 0 :
       if self.cholxyz.mean(axis=0)[2] < self.protxyz.mean(axis=0)[2] :
         self.cholleaf = "low"
       else :
         self.cholleaf = "upp"

    # Optionally, read a file of Lennard-Jones sigma parameters
    if sigmafile is None :
      self.sigmas = np.ones(self.protxyz.shape[0])*2.3
    else :
      self.sigmas = []
      with open(sigmafile,'r') as f :
        line = f.readline()
        while line :
          self.sigmas.append(float(line.strip()))
          line = f.readline()
      self.sigmas = np.array(self.sigmas)

  def plot(self,axis, leaflet, centroid=None, reverseY=True,sigmas=None,sidechain=False,
            colorscheme=None, drawchol=True, specialres=None, rotate2D=None) :
    """
    Plot the structure on a 2D grid

    Parameters
    ----------
    axis : Axis object
      the matplotlib Axis object to draw on
    leaflet : string
      the leaflet to draw
    reverseY : boolean, optional
      if the y coordinate should be reversed
      used to distinguish between an look-up or look-down view
    sigmas : numpy array, optional
      sizes of atoms
    sidechain : boolean, optional
      if we should draw sidechains
    colorscheme : list, optional
      colors to use if not pre-defined
    drawchol : boolean, optional
      if to draw cholesterols if they exists
    specialres : list of int
      residues to plot in black
    """
    if sigmas is None and self.sigmas is not None :
      sigmas = self.sigmas

    if specialres is not None:
      if colorscheme is None:
        colorscheme = -np.ones([self.pdbfile.xyz.shape[0],3])
      for res in specialres :
        resid = self.template.residues.index(res)
        for atom in self.pdbfile.residues[resid].atoms:
          colorscheme[atom.idx,:] = [0.05,0.05,0.05]

    grid = GridTemplate(np.array([[0.0,0.0,0.0],self.box]),resolution=0.5)
    centx = int(float(grid.size[0])/2.0)#+(35.0-centroid[0])*grid.resolution
    centy = int(float(grid.size[1])/2.0)#+(35.0-centroid[1])*grid.resolution
    #print centx*grid.resolution,centy*grid.resolution,centroid,-(35.0-centroid[1]),-(35.0-centroid[0])

    coords = np.array(self.protxyz,copy=True)
    if rotate2D is not None :
      m = self.protxyz.mean(axis=0)[:2]
      coords[:,:2] = rotate2D(coords[:,:2]-m)+m

    hpoints = np.zeros([len(self.helices),2])
    for i,h in enumerate(self.helices) :
      hlxxyz = coords[h[0]:h[1]+1,:]
      aid = np.arange(h[0],h[1]+1)
      if leaflet[0:3] == "upp" :
        sel = hlxxyz[:,2] >= self.box[2]/2
      else :
        sel = hlxxyz[:,2] <= self.box[2]/2
      if sel.sum() > 0 :
        idx = grid.indices(hlxxyz[sel,:])
        hpoints[i,0] = (idx[0].mean()-centx)*grid.resolution
        hpoints[i,1] = ((-idx[1]).mean()+centy)*grid.resolution

        hlxsig = sigmas[h[0]:h[1]+1]
        hlxflag = self.protflag[h[0]:h[1]+1]
        for ai,coord,sig,flag in zip(aid[sel],hlxxyz[sel,:],hlxsig[sel],hlxflag[sel]) :
          if flag == -1 and not sidechain : continue
          fac = -1
          if reverseY : fac = 1
          thiscolor = colors.color(i)
          if i == 0 : thiscolor = colors.color(11)
          thisalpha = 0.8
          if colorscheme is not None :
            if colorscheme[ai,0] > -1  :
              thiscolor = colorscheme[ai,:]
              thisalpha = 1.0
          cir = plt.patches.Circle((fac*(-coord[1]+centy*grid.resolution),coord[0]-centx*grid.resolution),
                radius=sig,fc=thiscolor,ec=colors.darker(thiscolor),alpha=thisalpha)
          axis.add_patch(cir)

    # Add annotation for each helix
    for i,pnt in enumerate(hpoints) :
      fac = -1
      if reverseY : fac = 1 #colors.color(i)
      axis.text(fac*pnt[1],pnt[0],"%d"%(i+1),bbox=dict(facecolor='white',alpha=0.9),color='k',fontsize=10)

    # Write out cholesterols
    if drawchol and leaflet[0:3] == self.cholleaf and self.cholxyz is not None and self.cholxyz.shape[0] > 0 :
      fac = -1
      if reverseY : fac = 1
      coords = np.array(self.cholxyz,copy=True)
      if rotate2D is not None :
        m = self.protxyz.mean(axis=0)[:2]
        coords[:,:2] = rotate2D(coords[:,:2]-m)+m
      axis.plot(fac*(-coords[::2,1]+centy*grid.resolution),coords[::2,0]-centx*grid.resolution,'xk',markeredgewidth=1.0)

def draw_colormap(figure,image, text=r'$-RT\ \ln(\rho/\rho_0)$', unittxt="$\mathrm{[kJ/mol]}$") :
  """
  Draw a colorbar associated with an image
  """
  cax = figure.add_axes([ 0.08,0.0, 1,1])
  if unittxt != "":
      cax.text(0.99,0.80, text)
      cax.text(1.03,0.77, unittxt)
  else:
      cax.text(0.99,0.77,text)
  hide_axis(cax)
  figure.colorbar(image,orientation='vertical',ax=cax,shrink=0.5,aspect=50)
  return cax

def draw_joint2d(axis,residues0,residues,contacts,logit=False) :
  """
  Draw a residue-residue contact joint probability plot as a contour

  Parameters
  ----------
  axis : Axis object
    the axis to draw on
  residues0 : numpy array
    the tick positions of the residues
  residues : numpy array
    the x-ray number of the residues
  contacts : numpy array
    the contact probabilities
  logit : boolean, optional
    if to plot log probabilities
  """

  sel = contacts>0
  if logit :
    logcont = np.zeros(contacts.shape)
    logcont[sel] = np.log(contacts[sel])
    C = axis.imshow(logcont,cmap=plt.cm.BuGn,origin="lower",interpolation="none")
  else :
    C = axis.imshow(contacts,cmap=plt.cm.BuGn,origin="lower",interpolation="none")
  excluded_im = np.ones([contacts.shape[0],contacts.shape[1],4])
  excluded_im[sel,3] = 0.0
  axis.imshow(excluded_im,origin="lower")

  sel = residues > 0
  axis.set_xticks(residues0[sel][::20])
  axis.set_yticks(residues0[sel][::20])
  axis.set_xticklabels(residues[sel][::20])
  axis.set_yticklabels(residues[sel][::20])
  axis.set_xlabel("Residue")
  axis.set_ylabel("Residue")

  return C

def hide_axis(axis) :
  """
  Hide the axes and frame of matplotlib axis object
  """
  axis.get_xaxis().set_visible(False)
  axis.get_yaxis().set_visible(False)
  axis.patch.set_alpha(0)
  axis.set_frame_on(False)

def plot_density_xray(axis,density,mat,minval,maxval,xray,leaflet,label,
    number=None, plotn=True, drawchol=True, specialres=None) :
  """
  Utility function to plot both density and x-ray structure

  Parameters
  ----------
  axis : Axis object
    the axis to plot on
  density : Density
    the density to plot
  mat : string
    the matrix of the Density to plot
  minval : float
    the minimum value of the image
  maxval : float
    the maximum value of the image
  xray : XrayDensity
    the xray structure to plot
  leaflet : string
    the leaflet to draw
  label : string
    a label to draw inside the image
  number : string, optional
    a number string to draw to identify this plot with, e.g. A), B)
  plotn : boolean, optional
    if to draw the average number of molecules on the density
  drawchol : boolean, optional
      if to draw cholesterols if they exists
  specialres : list of int
    residues to plot in black
  Returns
  -------
  Image :
    the Image of the density that was plotted, for adding colormap
  """

  if density == 0 :
    density = Density(np.zeros([70,70]),2)
    density.append(np.zeros([70,70]))
    density.average()
    mat = "average"
    minval = -1000
    maxval = 1000

  def rotate_mat(m):
    return np.rot90(m, 2)
  rotate2D =  None # if leaflet=="low" else rotate_mat
  c = density.centroid(rotate2D)

  if density is not None :
    im = density.plot(mat,axis,plt.cm.RdYlBu,minval,maxval=maxval ,reverseY=leaflet=="upp", rotate2D=rotate2D)
  else :
    im = None

  def rotate_coords(coords):
    angle = np.pi
    r = np.matrix( ((np.cos(angle),-np.sin(angle)), (np.sin(angle), np.cos(angle))) )
    return np.asarray(coords*r)
  rotate2D = None  #if leaflet=="low" else rotate_coords

  xray.plot(axis,leaflet, c, reverseY=leaflet=="upp",drawchol=drawchol, specialres=specialres, rotate2D=rotate2D)
  axis.set_xlim((-35,35))
  axis.set_xticks([-30,-20,-10,0,10,20,30])
  axis.set_yticks([-30,-20,-10,0,10,20,30])
  axis.set_xticklabels(["-3","-2","-1",0,"1","2","3"])
  axis.set_yticklabels(["-3","-2","-1",0,"1","2","3"])
  axis.set_ylim((-35,35))
  axis.text(-28,28,label,size=14)
  if number is not None : axis.text(-45,38,number)
  if plotn : axis.text(-28,-32,"<n>=%.1f"%(density.molcount.mean()),size=12)
  return im

################
# Misc routines
################

def load_template(mol) :
  """
  Load a ProteinTemplate object from a standard location
  """
  filename = os.path.join(PROT_INFO_PATH,"template_%s.txt"%mol)
  if os.path.isfile(filename) :
    return ProteinTemplate(filename)
  else :
    raise Exception("Invalid mol (%s) or file is missing (%s)"%(mol,filename))

def load_xray(mol, loadsigma=False, loadaa=False) :
  """
  Load an XrayDensity object from a standard location
  """
  template = load_template(mol)
  filename = os.path.join(PROT_INFO_PATH,"xray_%s.gro"%mol)
  if os.path.isfile(filename) :
    sigmafile = None
    # Optionally give the sigma file
    if loadsigma :
      sigmafile = os.path.join(PROT_INFO_PATH,"sigmas_%s"%mol)
      if not os.path.isfile(sigmafile) : sigmafile = None

    if loadaa :
      filename2 = os.path.join(PROT_INFO_PATH,"xray_%s-aa.gro"%mol)
      return XrayDensity(filename,template,sigmafile), pdb.PDBFile(filename2)
    else:
      return XrayDensity(filename,template,sigmafile)
  else :
    raise Exception("File is missing (%s)"%filename)

def read_rescontacts(folder, mol, percentile=90, expand=0, returnprobs=False):
  """
  Reads residue contacts from file and return a percentile
  """
  residues = []
  contactprob = []
  with open(os.path.join(folder,"%s_chol_6A_com_1d.txt"%mol), "r") as f :
    lines = f.readlines()
  for line in lines:
    res, p, d1, d2 = line.strip().split()
    contactprob.append(float(p))
    residues.append(int(res[3:]))
    for i in range(expand):
        residues.append(int(res[3:])-i)
        residues.append(int(res[3:])+1)
  residues = np.asarray(residues)
  contactprob = np.asarray(contactprob)
  sel = contactprob >= np.percentile(contactprob, percentile)
  if returnprobs:
    return residues[sel],contactprob[sel]
  else:
    return residues[sel]

def logical_expr(expr,*args) :
  """
  Simple parser of logical expression, to apply logical operators (not, and, or) on
  numpy arrays.

  Supported characters:
  1-9 : index into args list of numpy array, start from 1
  n : logical not operator
  a : logical and operator
  o : logical or operator

  Does not support parenthesis!
  """
  trace = ""
  for i,c in enumerate(expr) :
    if c not in ['n','a','o'] :
      trace = trace+c
      if i == len(expr) - 1 :
        trace = int(trace)
        return args[trace-1]
    elif c == "n" :
      return np.logical_not(logical_expr(expr[i+1:],*args))
    elif c == "a" :
      trace = int(trace)
      return np.logical_and(args[trace-1],logical_expr(expr[i+1:],*args))
    elif c == "o" :
      trace = int(trace)
      return np.logical_or(args[trace-1],logical_expr(expr[i+1:],*args))

def read_booleans(fileobj) :
  """
  Read a list of booleans from a file in binary format
  Returns False is there is nothing more to read from the file
  """
  line = fileobj.readline().strip()
  if line :
    return True,np.fromstring(line,dtype=np.uint8)
  else :
    return False,[]

def read_statefile(filename,every=1,block=None) :
  """
  Read a state file from disc
  This is a file where on each line there is a bit string
  """
  state = []
  with open(filename,'rb') as f :
    flag,list = read_booleans(f)
    while flag :
      state.append(list)
      flag,list = read_booleans(f)
  if every == 1 :
    state =  np.array(state,dtype=np.uint8)
  else :
    state = np.array(state,dtype=np.uint8)[::every]
  if block is not None:
      n = int(float(state.shape[0]) / float(block[0]) * block[1])
      return state[:n+1, :]
  else:
    return state

def write_booleans(fileobj,list) :
  """
  Write list of booleans to a file in binary format
  """
  bstr = struct.pack("?"*len(list),*list)
  #bstr = "".join("%d"%b for b in list)
  fileobj.write("%s\n"%bstr)
