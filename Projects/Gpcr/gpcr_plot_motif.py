# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to draw CRAC and CCM motifs on x-ray structures

No arguments are necessary, all structures are taken from standard locations
"""

import os
import sys

import numpy as np

import matplotlib
if not "DISPLAY" in os.environ or os.environ["DISPLAY"] == "" :
  matplotlib.use('Agg')
import matplotlib.pyplot as plt

import gpcr_lib

def _make_cra_colors(struct,colors=None) :

  if colors is None :
    colors = -np.ones([struct.xyz.shape[0],3])
  for i,res in enumerate(struct.residues) :
    if res.resname.strip() == "TYR" :
      left = -1
      for j in range(max(i-6,0),max(i-2,0)+1) :
        if struct.residues[j].resname.strip() in ["LEU","VAL"] :
          left = j
          break
      right = -1
      if left > -1 :
        for j in range(min(i+2,len(struct.residues)),min(i+6,len(struct.residues))+1) :
          if struct.residues[j].resname.strip() in ["LYS","ARG"] :
            right = j
            break
      if left > -1 and right > -1 :
        print "CRA: "+" ".join(res2.resname for res2 in struct.residues[left:right+1])
        for res2 in struct.residues[left:right+1] :
          for atom in res2.atoms :
            colors[atom.idx,:] = np.array([201.0/255.0,148.0/255.0,199.0/255.0])
  return colors

def _make_ccm_colors(struct,colors=None) :

  if colors is None :
    colors = -np.ones([struct.xyz.shape[0],3])
  for i,res in enumerate(struct.residues) :
    if res.resname.strip() in ["TYR","TRP"] :
      left = -1
      for j in range(max(i-11,0),max(i-8,0)+1) :
        if struct.residues[j].resname.strip() in ["LYS","ARG"] :
          left = j
          break
      right = -1
      if left > -1 :
        for j in range(max(i-4,0),max(i-4,0)+1)  :
          if struct.residues[j].resname.strip() in ["VAL","LEU","ILE"] :
            right = j
            break
      if left > -1 and right > -1 :
        print "CCM: "+" ".join(res2.resname for res2 in struct.residues[left:i+1])
        for res2 in struct.residues[left:i+1] :
          for atom in res2.atoms :
            colors[atom.idx,:] = np.array([216.0/255.0,179.0/255.0,101.0/255.0])
  return colors

if __name__ == '__main__' :

  respath = sys.argv[1]
  mols = "b2 b2_a a2a a2a_a".split()
  numbers = "A) B) C) D) E) F) G) H)".split()
  fig = plt.figure(1,figsize=(8,12))

  for i, mol in enumerate(mols) :

    xray = gpcr_lib.load_xray(mol)
    pdb = xray.pdbfile

    # Make colour pattern due to pattern
    colors = _make_cra_colors(pdb)
    colors = _make_ccm_colors(pdb,colors)
    restocolor = gpcr_lib.read_rescontacts(respath, mol)

    a = fig.add_subplot(len(mols),2,i*2+1,aspect="equal")
    xray.plot(a,"low",reverseY=False,sidechain=True,colorscheme=colors,drawchol=False, specialres=restocolor)
    a.text(-35,33,numbers[(i*2)])
    a.text(-25,25,"Intra.",size=14)
    a.set_xticklabels([])
    a.set_yticklabels([])
    a.set_xlim((-30,30))
    a.set_ylim((-30,30))
    if i == 0 :
      c = plt.Circle([-25,-27],radius=2.0,fc=[201.0/255.0,148.0/255.0,199.0/255.0],ec=[201.0/255.0,148.0/255.0,199.0/255.0])
      a.add_patch(c)
      a.text(-22,-28.5,"CRAC",size=14)
      c = plt.Circle([-5,-27],radius=2.0,fc=[216.0/255.0,179.0/255.0,101.0/255.0],ec=[216.0/255.0,179.0/255.0,101.0/255.0])
      a.add_patch(c)
      a.text(-1,-28.5,"CCM",size=14)
      c = plt.Circle([15,-27],radius=2.0,fc='k',ec='k')
      a.add_patch(c)
      a.text(20,-28.5,"Res.",size=14)

    a = fig.add_subplot(len(mols),2,i*2+2,aspect="equal")
    xray.plot(a,"upp",reverseY=True,sidechain=True,colorscheme=colors,drawchol=False)
    a.text(-35,33,numbers[(i*2)+1])
    a.text(-25,25,"Extra.",size=14)
    a.set_xticklabels([])
    a.set_yticklabels([])
    a.set_xlim((-30,30))
    a.set_ylim((-30,30))

  fig.savefig("motifs.png",format="png")
