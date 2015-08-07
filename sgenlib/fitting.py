# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Routines to fit structures
"""

import math
import sys

import numpy as np


def dofit(ref,mob,xyz) :
  """
  Performing fitting of mobile coordinates on top of references
  coordinates and then move a full set of coordinates
  
  Parameters
  ----------
  ref : numpy.ndarray
    the reference coordinates
  mob : numpy.ndarray
    the mobile coordinates
  xyz : numpy.ndarray
    the full set of coordinates to move

  Returns
  -------
  numpy.ndarray
    the new coordinates
  """
  ref_com = np.average(ref,axis=0)
  mob_com = np.average(mob,axis=0)
  
  ref = ref - ref_com
  mob = mob - mob_com
  xyz = xyz - mob_com

  rotmat = kabschFit(ref,mob)
  xyz = rotate(xyz,rotmat)

  xyz = xyz + ref_com

  return xyz


def kabschFit(ref,mob) :
  """
  Performs a fit of coordinate using the Kabsch algorithm
  
  Parameters
  ----------
  ref : numpy.ndarray
    the reference coordinates
  mob : numpy.ndarray
    the mobile coordinates

  Returns
  -------
  numpy.ndarray
    the rotation matrix
  """
  # Form R matrix
  r = np.zeros((3,3))
  for i in range(3) :
    for j in range(3) :
      r[i,j] = np.sum(np.multiply(ref[:,i],mob[:,j]))
  
  # Form Rtr matrix
  rtr = np.zeros((3,3))
  for i in range(3) :
    for j in range(3) :
      rtr[i,j] = np.sum(np.multiply(r[i,:],r[j,:]))

  # Find eigenvalues and eigenvectors or rtr
  (eigenval,eigenvec) = np.linalg.eig(rtr)
  idxs = eigenval.argsort()
  eigenval = eigenval[idxs]
  eigenvec = eigenvec[:,idxs]

  eigenvec[:,0] = np.cross(eigenvec[:,1],eigenvec[:,2])

  # Contruct the B matrix  
  b = np.zeros((3,3))
  for i in range(3) :
    for j in range(3) :
      b[i,j]  = np.sum(np.multiply(r[:,i],eigenvec[:,j])/np.sqrt(np.abs(eigenval[j])))
  b[:,1] = b[:,1] / np.linalg.norm(b[:,1])
  b[:,2] = b[:,2] / np.linalg.norm(b[:,2])
  b[:,0] = np.cross(b[:,1],b[:,2])

  # Contruct the rotation matrix
  rotmat = np.zeros((3,3))
  for i in range(3) :
    for j in range(3) :
      rotmat[i,j] = np.sum(np.multiply(eigenvec[i,:],b[j,:]))

  return rotmat

def rotate(xyz,rotmat) :
  """
  Rotate a set of coordinates
  
  Parameters
  ----------
  xyz : numpy.ndarray
    coordinates to rotate
  rotmat : numpy.ndarray
    rotation matrix
    
  Returns
  -------
  numpy.ndarray
    the rotate coordinates
  """
  xyz2 = np.zeros(xyz.shape)
  for i in range(xyz.shape[0]) :
    for j in range(3) :
      for k in range(3) :
        xyz2[i,j] = xyz2[i,j]+xyz[i,k]*rotmat[j,k]
  return xyz2

"""
#
# Unclear if the below functions works as expected, debugging needed
#

def make_rotmat(angle,vnorm,coord,transpose):
    ### From: http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/ Section 6.2
    u,v,w = vnorm[0],vnorm[1],vnorm[2]
    a,b,c = coord[0],coord[1],coord[2]
    makerotmat = np.zeros([4,4])
    makerotmat[0,0] = (math.pow(u,2)+(math.pow(v,2)+math.pow(w,2))*math.cos(angle))
    makerotmat[0,1] = (u*v*(1-math.cos(angle))-w*math.sin(angle))
    makerotmat[0,2] = (u*w*(1-math.cos(angle))+v*math.sin(angle))
    makerotmat[0,3] = (a*(math.pow(v,2)+math.pow(w,2))-u*(b*v+c*w))*(1-math.cos(angle))+(b*w-c*v)*math.sin(angle)
    makerotmat[1,0] = (u*v*(1-math.cos(angle))+w*math.sin(angle))
    makerotmat[1,1] = (math.pow(v,2)+(math.pow(u,2)+math.pow(w,2))*math.cos(angle))
    makerotmat[1,2] = (v*w*(1-math.cos(angle))-u*math.sin(angle))
    makerotmat[1,3] = (b*(math.pow(u,2)+math.pow(w,2))-v*(a*u+c*w))*(1-math.cos(angle))+(c*u-a*w)*math.sin(angle)
    makerotmat[2,0] = (u*w*(1-math.cos(angle))-v*math.sin(angle))
    makerotmat[2,1] = (v*w*(1-math.cos(angle))+u*math.sin(angle))
    makerotmat[2,2] = (math.pow(w,2)+(math.pow(u,2)+math.pow(v,2))*math.cos(angle))
    makerotmat[2,3] = (c*(math.pow(u,2)+math.pow(v,2))-w*(a*u+b*v))*(1-math.cos(angle))+(a*v-b*u)*math.sin(angle)
    makerotmat[3,0] = 0.0
    makerotmat[3,1] = 0.0
    makerotmat[3,2] = 0.0
    makerotmat[3,3] = 1.0
    if transpose :
      return makerotmat.transpose()
    else :
      return makerotmat

def rotateAround(atom1,atom2,angle,xyz,transpose) :

  vec = atom2 - atom1
  vlen = np.sqrt((vec**2).sum())
  vec = vec / vlen
  rotmat = make_rotmat(angle,vec,atom2,transpose)
  xyz = xyz - atom2
  xyz = rotate(xyz,rotmat)
  xyz = xyz + atom2
  return xyz
"""
    
