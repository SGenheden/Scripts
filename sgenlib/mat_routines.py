# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Helper routines for density processing, from the Gpcr project originally
"""

import numpy as np

def read_simple(filename,lownam,uppnam) :

  mat = np.load(filename)
  return mat[lownam],mat[uppnam]

def read_and_scale(filename,lownam,uppnam,scaling=None) :

  lowmat,uppmat = read_simple(filename,lownam,uppnam)

  if scaling["flag"] != None :
    lowsel = lowmat!=0.000
    uppsel = uppmat!=0.000
    if scaling["flag"] > 0 :
      uppmat[uppsel] = 1.987*0.3*np.log(uppmat[uppsel]/scaling["factor"][1])
      lowmat[lowsel] = 1.987*0.3*np.log(lowmat[lowsel]/scaling["factor"][0])
    else :
      uppmat[uppsel] = uppmat[uppsel]/scaling["factor"][1]
      lowmat[lowsel] = lowmat[lowsel]/scaling["factor"][0]

  return lowmat,uppmat

def read_and_filter(filename,lownam,uppnam,lowfilter,uppfilter,scaling=None) :

  lowmat,uppmat=read_and_scale(filename,lownam,uppnam,scaling)

  if lowfilter != None and uppfilter != None :
    lowmat[lowfilter] = 0.0
    uppmat[uppfilter] = 0.0

  return lowmat,uppmat

def mol_dist(filename,nam) :


  MAT_SIZE = 70.0*70.0
  mat = np.load(filename)
  return mat[nam][:,0].mean()/float(MAT_SIZE),mat[nam][:,1].mean()/float(MAT_SIZE)

def set_filter(filter_opt,scaling,index) :

  if filter_opt["density"] == None :
    return None,None

  lowmat,uppmat = read_and_scale(filter_opt["density"][index],filter_opt["lower"],filter_opt["upper"],scaling)

  return lowmat < filter_opt["threshold"],uppmat < filter_opt["threshold"]
