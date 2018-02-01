# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Routines to geometrical calculations
"""

import math

import numpy as np

def sphere_surf_rand():

  u1 = 100
  u2 = 1000
  while u1*u1 + u2*u2 >= 1 :
    u1 = np.random.uniform(low=-1,high=1)
    u2 = np.random.uniform(low=-1,high=1)
  x = 2*u1*np.sqrt(1-u1*u1-u2*u2)
  y = 2*u2*np.sqrt(1-u1*u1-u2*u2)
  z = 1 - 2*(u1*u1 + u2*u2)
  return np.asarray([x,y,z])

def sphere_rand(r) :
  """
  Random number on sphere with radius r
  """
  z = np.random.uniform(low=-r,high=r)
  phi = np.random.uniform(low=0,high=2.0*np.pi)
  theta = np.arcsin(z/r)
  x=r*np.cos(theta)*np.cos(phi)
  y=r*np.cos(theta)*np.sin(phi)
  return np.array([x,y,z])

def vnorm(v) :
  """
  Returns a normalised vector
  """
  vlen = np.sqrt((v**2).sum())
  return v / vlen

def build_xyz(bndatm,angatm,dihedatm,blen,ang,dih) :
  """
  Build Cartesian coordinates from internal coordinates
  """

  ang = ang*np.pi/180.0
  dih = dih*np.pi/180.0
  blen = blen

  #print bndatm[0],bndatm[1],bndatm[2],angatm[0],angatm[1],angatm[2],dihedatm[0],dihedatm[1],dihedatm[2],dih,ang,blen

  d1 = dihedatm - bndatm
  d2 = angatm   - bndatm

  vy = vnorm(np.cross(d1,d2))
  vx = vnorm(np.cross(d2,vy))
  vz = vnorm(np.cross(vx,vy))

  xbs= blen * np.sin(ang) * np.cos(dih)
  ybs = blen * np.sin(ang) * np.sin(dih)
  zbs = -blen * np.cos(ang)

  xyz = np.zeros(3)
  xyz[0] = vx[0]*xbs + vy[0]*ybs + vz[0]*zbs + bndatm[0]
  xyz[1] = vx[1]*xbs + vy[1]*ybs + vz[1]*zbs + bndatm[1]
  xyz[2] = vx[2]*xbs + vy[2]*ybs + vz[2]*zbs + bndatm[2]

  return xyz

def angle(v1,v2) :
  """
  Calculates the angle between two vectors
  """
  a = (v1*v2).sum()/(np.sqrt((v1**2).sum())*np.sqrt((v2**2).sum()))
  if a > 1.0 :
    a = 0.0
  elif a < -1 :
    a = np.pi
  else :
    a = np.arccos(a)
  return a

def angle_atms(a1,a2,a3) :
  """
  Calculates the angle between three vectors/atoms
  """
  v1 = a2 - a1
  v2 = a2 - a3
  return angle(v1,v2)

def dihedral(a1,a2,a3,a4) :
  """
  Calculates the dihedral angle between four vectors/atoms
  """
  v1 = a2 - a1
  v2 = a3 - a2
  v3 = a4 - a3
  n1 = np.cross(v1,v2)
  n2 = np.cross(v2,v3)
  if np.dot(v3, np.cross(v1, v2)) <= 0.0 :
    return -angle(n1,n2)
  else :
    return angle(n1,n2)


def dihedral_protoms(a1,a2,a3,a4) :
  """
  Calculates a dihedral angles exactly as in ProtoMS
  """
  v21 = a1 - a2
  v23 = a3 - a2
  v32 = -v23
  v34 = a4 - a3

  n1 = vnorm(np.cross(v21,v23))
  n2 = vnorm(np.cross(v32,v34))

  phi = np.dot(n1,n2)

  if phi > 1.0 :
    phi = 1
  elif phi < -1.0 :
    phi = -1

  phi = np.arccos(phi)
  n3 = vnorm(np.cross(n1,v23))
  ang = np.dot(n2,n3)
  if ang > 1.0 :
    ang = 1.0
  elif ang < -1 :
    ang = -1.0
  ang = np.arccos(ang)

  if ang < np.pi/2.0 : phi = 2*np.pi-phi

  return phi

def multi2RB(forces,multiplicities,phases) :
  """
  Converts a period dihedral angle to Ryckart-Bell dihedral
  """
  V = 5*[0.0]
  C = 5*[0.0]

  for i in range(len(forces)) :
    force = forces[i]
    period = multiplicities[i]
    phase = phases[i]

    if force > 0 : V[period] = 2*force
    if period == 1 :
      C[0] += 0.5 * V[period]
      if phase == 0:
        C[1] -= 0.5 * V[period]
      else:
        C[1] += 0.5 * V[period]
    elif period == 2:
      if phase == 180:
        C[0] += V[period]
        C[2] -= V[period]
      else:
        C[2] += V[period]
    elif period == 3:
      C[0] += 0.5 * V[period]
      if phase == 0:
        C[1] += 1.5 * V[period]
        C[3] -= 2 * V[period]
      else:
        C[1] -= 1.5 * V[period]
        C[3] += 2 * V[period]
    elif period == 4:
      if phase == 180:
        C[2] += 4 * V[period]
        C[4] -= 4 * V[period]
      else:
        C[0] += V[period]
        C[2] -= 4 * V[period]
        C[4] += 4 * V[period]

  return C


def sphericity(xyz) :

    xyz2 = xyz - xyz.mean(axis=0)
    tensor = np.zeros((3, 3))
    for x in range(xyz2.shape[0]):
        tensor += np.outer(xyz2[x, :], xyz2[x, :])
    eigval = np.linalg.eigvalsh(tensor)
    shape = (3.0 / 2.0) * 3 * np.var(eigval) / np.sum(eigval)**2
    return eigval, shape

def centerOfMass(xyz,masses) :
  """
  Calculates the centre of mass of a set of coordinates
  """
  return np.sum(xyz*masses[:,np.newaxis],axis=0)/np.sum(masses)

def moment_of_inertia(xyz,masses):
  """
  Calculates the moment of inertia of a set of coordinates
  """
  # Convert to local coordinates
  recenteredpos = xyz - centerOfMass(xyz,masses)
  values = zip(masses, recenteredpos)
  # Create the inertia tensor
  # m_i = mass of atom i
  # (x_i, y_i, z_i) = pos of atom i
  # Ixx = sum(m_i*(y_i^2+z_i^2)); Iyy = sum(m_i*(x_i^2+z_i^2)); Izz = sum(m_i*(x_i^2+y_i^2))
  # Ixy = Iyx = -1*sum(m_i*x_i*y_i)
  # Ixz = Izx = -1*sum(m_i*x_i*z_i)
  # Iyz = Izy = -1*sum(m_i*y_i*z_i)
  Ixx = reduce(lambda t,a: t+a[0]*(a[1][1]*a[1][1]+a[1][2]*a[1][2]), values, 0.)
  Iyy = reduce(lambda t,a: t+a[0]*(a[1][0]*a[1][0]+a[1][2]*a[1][2]), values, 0.)
  Izz = reduce(lambda t,a: t+a[0]*(a[1][0]*a[1][0]+a[1][1]*a[1][1]), values, 0.)
  Ixy = Iyx = -1*reduce(lambda t,a: t+a[0]*a[1][0]*a[1][1], values, 0.)
  Ixz = Izx = -1*reduce(lambda t,a: t+a[0]*a[1][0]*a[1][2], values, 0.)
  Iyz = Izy = -1*reduce(lambda t,a: t+a[0]*a[1][1]*a[1][2], values, 0.)
  return np.array([[Ixx, Ixy, Ixz],[Iyx, Iyy, Iyz],[Izx, Izy, Izz]])

def principal_axes(moi):
  """
  Calculates the principal axis of a moment of inertia vector
  """
  from numpy.linalg import eig
  eigenval, eigenvec = eig(moi)
  indices = np.argsort(eigenval)
  return eigenvec[:,indices].T

def rotate_coords(coords,angles=None):
    """
    Rotate coordinates using Euler angles
    """
    if angles is None :
        angles = [np.random.uniform(-np.pi,np.pi),
                    np.random.uniform(0,np.pi),
                    np.random.uniform(-np.pi,np.pi)]

    # pich - roll - yawn convention
    X =  np.mat([[1.0, 0.0, 0.0],
                    [0.0, np.cos(angles[0]), np.sin(angles[0])],
                    [0.0, -np.sin(angles[0]), np.cos(angles[0])]])
    Y = np.mat([[np.cos(angles[1]), 0.0, -np.sin(angles[1])],
                    [0.0, 1.0, 0.0],
                    [np.sin(angles[1]), 0.0, np.cos(angles[1])]])
    Z = np.mat([[np.cos(angles[2]), np.sin(angles[2]),0.0],
                    [-np.sin(angles[2]), np.cos(angles[2]),0.0],
                    [0.0, 0.0, 1.0]])

    newmat = np.mat(coords-coords.mean(axis=0)).T
    rotated = Z*Y*X*newmat
    return np.asarray(rotated.T + np.mean(coords,axis=0))

def rotaxis(a,b):
  """
  Rotates a vector around an axis
  """
  c = np.cross(a,b)
  return vnorm(c)

def rotation_matrix(angle, direction, point=None):
  """
  Constructs a rotation matrix
  """
  sina = math.sin(angle)
  cosa = math.cos(angle)
  direction = vnorm(direction[:3])
  # rotation matrix around unit vector
  R = np.array(((cosa, 0.0,  0.0),
                   (0.0,  cosa, 0.0),
                   (0.0,  0.0,  cosa)), dtype=np.float64)
  R += np.outer(direction, direction) * (1.0 - cosa)
  direction *= sina
  R += np.array((( 0.0,         -direction[2],  direction[1]),
                    ( direction[2], 0.0,          -direction[0]),
                    (-direction[1], direction[0],  0.0)),
                   dtype=np.float64)
  M = np.identity(4)
  M[:3, :3] = R
  if point is not None:
      # rotation not around origin
      point = numpy.array(point[:3], dtype=np.float64, copy=False)
      M[:3, 3] = point - np.dot(R, point)
  return M

def vlen3(vec) :
  return np.sqrt((vec**2).sum())

def majorminor_axis(xyz) :

  def posneg_vec(vec1,vec2,posvec,negvec,npos,nneg,vlen_flag) :

    vlen = vlen3(vec1)*vlen3(vec2)
    vec = (vec1*vec2).sum()
    #if vlen_flag == 1 :
    #  vec = vec / vlen
    ang = np.arccos(vec/vlen)*180.0/np.pi
   # if vlen_flag == 0 :
   #   ang = ang / vlen
    if ang > 360.0 : ang = ang - 360.0
    if ang < 0.0 : ang = ang + 360.0
    if ang >= -90 and ang < 90 :
      posvec = posvec + vec2
      npos = npos + 1 #3
    else :
      negvec = negvec + vec2
      nneg = nneg + 1 # 3
    return posvec,negvec,npos,nneg

  def assign_vec(posvec,negvec,npos,nneg) :

    vec = np.zeros(3)
    if npos > 0 : posvec = posvec / float(npos)
    if nneg > 0 : negvec = negvec / float(nneg)

    poslen2 = (posvec**2).sum()
    neglen2 = (negvec**2).sum()

    if poslen2 > neglen2 :
      if poslen2 != 0.0 :
        vec = vnorm(posvec)
      else :
        vec[0] = 1.0
    else :
      if neglen2 != 0.0 :
        vec = vnorm(negvec)
      else :
        vec[0] = 1.0
    return vec

  center = xyz.mean(axis=0)
  major = np.zeros(3)
  minor_tmp = np.zeros(3)

  # Major axis
  posvec = np.zeros(3)
  negvec = np.zeros(3)
  npos = 0
  nneg = 0
  for i in range(xyz.shape[0]) :
    vec = xyz[i,:] - center
    posvec,negvec,npos,nneg = posneg_vec(np.array([1.0,0.0,0.0]),vec,posvec,negvec,npos,nneg,1)

  major = assign_vec(posvec,negvec,npos,nneg)

  # Minor axis
  posvec = np.zeros(3)
  negvec = np.zeros(3)
  npos = 0
  nneg = 0
  for i in range(xyz.shape[0]) :
    vec = xyz[i,:] - center
    posvec,negvec,npos,nneg = posneg_vec(major,vec,posvec,negvec,npos,nneg,0)

  minor_tmp = assign_vec(posvec,negvec,npos,nneg)

  majlen2 = (major**2).sum()
  minlen2 = (minor_tmp**2).sum()

  if majlen2 < 1E-10 :
    major = np.array([1.0,0.0,0.0])
  elif minlen2 < 1E-10 :
    minor_tmp = np.array([0.0,1.0,0.0])

  if minlen2 > majlen2 :
    temp = np.array(major,copy=True)
    major = np.array(minor_tmp,copy=True)
    minor_tmp = np.array(temp,copy=True)

  if np.abs((major*minor_tmp).sum()) > 0.9 :
    major_swap = np.zeros(3)
    major_swap[0] = major[0]
    major_swap[1] = major[2]
    major_swap[2] = major[1]
    minor = vnorm(np.cross(major,major_swap))
  else :
    minor = vnorm(np.array(minor_tmp,copy=True))
    alpha = (major*minor).sum()
    beta = ((minor**2).sum() - (alpha**2))
    beta = np.sqrt(beta)
    minor = (minor - (alpha*major)) / beta
    minor = vnorm(minor)

  return major,minor
