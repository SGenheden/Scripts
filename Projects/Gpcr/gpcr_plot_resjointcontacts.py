# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to plot residue joint contact probability

Examples
--------
gpcr_plot_rescontacts.py -f r1_md3_en_fit_joint.npz -m ohburr -l oh --mol b2
"""

import os
import argparse
import sys

import numpy as np
import matplotlib
if not "DISPLAY" in os.environ or os.environ["DISPLAY"] == "" :
  matplotlib.use('Agg')
import matplotlib.pyplot as plt

import pycontacts
import gpcr_lib
from sgenlib import colors
from sgenlib import pdb 

def _resjointcontacts(filename,label,mat,repeats,out,mol) :
  """
  Main analysis routine
  
  Parameters
  ----------
  filename : string
    file to analyse
  label : string
    label for the group
  mat : string
    matrix identifier
  repeats : list of string
    replacement pattern for multiple repeats
  out : string
    output prefix
  mol : string
    protein identifier
  """
  
  # Load the protein template to obtain residue information
  template = gpcr_lib.load_template(mol)
  residues = np.array(template.residues)
  residues0 = np.arange(1,residues.shape[0]+1)
  codes = np.array(template.codes)
  names = np.array(template.names)

  pcontacts = []
  npz = np.load(filename)
  pjoint0 = npz["joint"+mat]
  pjoint = np.zeros([len(repeats),pjoint0.shape[0],pjoint0.shape[0]])
  pjoint[0,:,:] = pjoint0 
        
  # Do the same for multiple repeats and average over them
  if repeats is not None :
    for ri,r in enumerate(repeats[1:],1) :
      filename2 = filename.replace(repeats[0],r)
      npz = np.load(filename2)
      pjoint[ri,:,:] = npz["joint"+mat]
      
    pjoint_std = pjoint.std(axis=0)/np.sqrt(len(repeats))
    f2d = plt.figure(2,tight_layout=True) 
    C = gpcr_lib.draw_joint2d(f2d.gca(),residues0,residues,pjoint_std)
    f2d.colorbar(C)
    f2d.savefig("%s_%s_2d_std.png"%(out,label),format="png")   
      
  # Draw a 2D residue-residue joint probability plot  
  f2d = plt.figure(1,tight_layout=True) 
  C = gpcr_lib.draw_joint2d(f2d.gca(),residues0,residues,pjoint.mean(axis=0))
  f2d.colorbar(C)
  f2d.gca().text(1.04,1.01,"p(A,B)",transform=f2d.gca().transAxes)
  f2d.savefig("%s_%s_2d.png"%(out,label),format="png") 
  
  f2d = plt.figure(3,tight_layout=True) 
  C = gpcr_lib.draw_joint2d(f2d.gca(),residues0,residues,pjoint.mean(axis=0),logit=True)
  f2d.colorbar(C)
  f2d.gca().text(1.02,1.02,"ln p(A,B)",transform=f2d.gca().transAxes)
  f2d.savefig("%s_%s_2d_log.png"%(out,label),format="png") 

if __name__ == '__main__' :

  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Plotting residue contacts")
  parser.add_argument('-f','--files',nargs='+',help="a list of input files.",default=[])
  parser.add_argument('-l','--labels',nargs='+',help="a label for each input file.",default=[])
  parser.add_argument('-m','--mat',nargs='+',help="the matrix to plot for each input file.",default=[])
  parser.add_argument('-o','--out',help="the output prefix.",default="rescontacts")
  parser.add_argument('--mol',choices=["b2","a2a","b2_a","a2a_a"],help="the protein molecules, should be either 'b2' or 'a2a'",default="b2")
  parser.add_argument('--repeats',nargs="+",help="replacement pattern for multiple repeats",default=["r1_","r2_","r3_","r4_","r5_"]) 
  args = parser.parse_args()
  
  for filename,label,mat in zip(args.files,args.labels,args.mat) :
    _resjointcontacts(filename,label,mat,args.repeats,args.out,args.mol)
