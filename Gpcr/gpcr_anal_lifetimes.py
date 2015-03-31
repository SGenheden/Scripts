# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to analyse state files in order to compute lifetimes of molecules or other
states

The state series can be files on disc or logical combinations of already open state
series. It will analyse group of files, each group with a number of repeats. 

Median and maximum lifetime will be written out standard output, as well as
average and standerd error in case of multiple repeats

Examples
--------
gpcr_anal_lifetimes.py -f r1_md3_en_fit_chol.mstate.6.dat r1_md3_en_fit_chol.mid.dat 1a2
                          r1_md3_en_fit_chol.resstate.6.dat
                       -l chol/mol chol/bur chol-on/bur chol/hlx --helical 4 --mol b2
                    
will read cholesterol in contact and buried cholesterol, and combine it to also analyse
buried cholesterols in contact. The fourth file exemplifies the helical lifetime analysis.
"""

import os
import sys
import argparse

import numpy as np

import pycontacts
import gpcr_lib

SEP = "\n"+" "*18


def _make_helical(state,mol) :
  """
  Convert a residual state matrix to a
  helical state matrix by looking at the residues
  in the helices.
  
  Parameters
  ----------
  state : nxR Numpy array
    the state matrix, 
    n is the number of snapshots
    R is the number of residues
  mol : string
    the identifier of the molecule

  Returns
  -------
  nxH Numpy array
    the helical state matrix, H is the number of helices 
  """
  helices = gpcr_lib.load_template(mol).rhelices
  hstate = np.zeros([state.shape[0],len(helices)],dtype=np.uint8)
  for i,h in enumerate(helices) :
    hstate[:,i] = np.any(state[:,h[0]-1:h[1]],axis=1)
  return hstate

def _anal_lifetimes(files,labels,helical,mol,time) :
  """
  Main work routine to analyse life times
  
  Parameters
  ----------
  files : list of string
    the state files to analyse
  labels : list of string
    a label for each state file
  helical : list of integers
    indicator to do helical transformation for some files
  mol : string
    identifier of molecule
  time : float
    the total simulation time
    
  Returns
  -------
  list of Numpy arrays
    the results 
  """
  
  # Read in each state file or do a logical transformation
  states = []
  for filename,label in zip(files,labels) :
    if os.path.isfile(filename) :
      states.append(gpcr_lib.read_statefile(filename))
    else :
      states.append(gpcr_lib.logical_expr(filename.replace(" ",""),*states))

  # Conversion factor from snapshot lifetime to ns lifetime
  ns_per_snapshot = time / float(states[0].shape[0])

  # Perform helical transformation     
  if helical is not None :
    for h in helical :
      states[h-1] = _make_helical(states[h-1],mol)   
     
  # Concatenate all state matrices and do lifetime analysis
  all_states = np.concatenate(states,axis=1)
  life_av,life_max = pycontacts.lifetime(all_states)   
        
  # Write out statistics
  print "%15s\t%8s\t%8s"%("","Median","Max")
  nused = 0
  results = []
  for i,(state,label) in enumerate(zip(states,labels),1) : 
    life_av_s = life_av[nused:nused+state.shape[1]]*ns_per_snapshot
    life_max_s = life_max[nused:nused+state.shape[1]]*ns_per_snapshot
    nused = nused+state.shape[1]
    if helical is not None and i in helical :
      print "%-18s%s"%(label,SEP.join("\t%8.3f\t%8.3f"%(a,m) for (a,m) in zip(life_av_s,life_max_s))) 
      results.append(np.concatenate([life_av_s,life_max_s]))
    else :
      print "%-18s\t%8.3f\t%8.3f"%(label,np.median(life_av_s),life_max_s.mean())
      results.append(np.array([np.median(life_av_s),life_max_s.mean()]))
      
  return results

if __name__ == '__main__' :

  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Analysing lifetimes from state files")
  parser.add_argument('-f','--files',nargs='+',help="a list of input files.",default=[])
  parser.add_argument('-l','--labels',nargs='+',help="a label for each state.",default=[])
  parser.add_argument('--helical',nargs='+',type=int,help="flag to perform helical transformation for a state file")
  parser.add_argument('--mol',choices=["b2","a2a","b2_a","a2a_a"],help="the protein molecules, should be either 'b2' or 'a2a'",default="b2")
  parser.add_argument('--time',type=float,help="total simulation time in ns",default=50000) 
  parser.add_argument('--repeats',nargs="+",help="replacement pattern for multiple repeats",default=["r1_","r2_","r3_","r4_","r5_"]) 
  args = parser.parse_args()
  
  res0 = _anal_lifetimes(args.files,args.labels,args.helical,args.mol,args.time)

  if args.repeats is not None :
    # Allocate arrays for total statistics
    results = []
    for r0 in res0 :
      results.append(np.zeros([r0.shape[0],len(args.repeats)]))
      results[-1][:,0] = r0
      
    for ri,r in enumerate(args.repeats[1:],1) :
      files2 = [f.replace(args.repeats[0],r) for f in args.files]
      print ""
      res = _anal_lifetimes(files2,args.labels,args.helical,args.mol,args.time)
      # Accumulate to the total statistics
      for i,r in enumerate(res) :
        results[i][:,ri] = r
        
    # Write out statistics over repeats    
    print "\n%15s\t%8s\t%8s"%("","Tot-Med.","Tot-Max")
    for i,(label,r) in enumerate(zip(args.labels,results),1) :
      rav = r.mean(axis=1) 
      rstd = r.std(axis=1)/np.sqrt(r.shape[1]) 
      if args.helical is not None and i in args.helical :
        l2 = rav.shape[0] / 2
        print "%-18s%s"%(label,SEP.join("\t%8.3f\t%8.3f\t%8.3f\t%8.3f"%(a,au,m,mu) for (a,m,au,mu) in zip(rav[:l2],rav[l2:],rstd[:l2],rstd[l2:]))) 
      else :
        print "%-18s\t%8.3f\t%8.3f\t%8.3f\t%8.3f"%(label,rav[0],rstd[0],rav[1],rstd[1])   