# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to analyse state files in order to count molecules 

The state series can be files on disc or logical combinations of already open state
series. It will analyse group of files, each group with a number of repeats. 

Average and standard deviation of the counts will be written out to standard output,
as well as average and standard error in case of multiple repeats

Examples
--------
gpcr_anal_counts.py -f r1_md3_en_fit_chol.leaflet.dat n1 
                    -l chol/intra. chol/extra. 
                    
will read cholesterols in the lower leaflet from disc and will negate this to also
analyse cholesterols in the upper leaflet.
"""

import os
import argparse

import numpy as np

import gpcr_lib

def _anal_contacts(filenames,labels) :
  """
  Main analysis routine
  
  Parameters
  ----------
  filenames : list of strings
    the files to read and analyse
  labels : list of strings
    the labels for the files

  Returns
  -------
  list of numpy array
    the produce result
  """
  
  # Read in state file from disc or combine them using expression evaluation
  states = []
  for filename,label in zip(filenames,labels) :
    if os.path.isfile(filename) :
      states.append(gpcr_lib.read_statefile(filename))
    else :
      states.append(gpcr_lib.logical_expr(filename.replace(" ",""),*states))
     
  # Count the number of on states and write out statistics      
  print "%18s\t%8s\t%8s"%("","Mean","Std")
  results = []
  for state,label in zip(states,labels) : 
    counton = state.sum(axis=1)
    print "%-18s\t%8.3f\t%8.3f"%(label,counton.mean(),counton.std())  
    results.append(np.array([counton.mean(),counton.std()]))
    
  return results

if __name__ == '__main__' :

  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Analysing counts from state files")
  parser.add_argument('-f','--files',nargs='+',help="a list of input files.",default=[])
  parser.add_argument('-l','--labels',nargs='+',help="a label for each state.",default=[])
  parser.add_argument('--repeats',nargs="+",help="replacement pattern for multiple repeats",default=["r1_","r2_","r3_","r4_","r5_"]) 
  args = parser.parse_args()
  
  res0 = _anal_contacts(args.files,args.labels)
  
  if args.repeats is not None :
    # Create arrays for average calculations
    results = []
    for r0 in res0 :
      results.append(np.zeros(len(args.repeats)))
      results[-1][0] = r0[0]
       
    for ri,r in enumerate(args.repeats[1:],1) :
      files2 = [f.replace(args.repeats[0],r) for f in args.files]
      print ""
      res = _anal_contacts(files2,args.labels)
      # Accumulate to average arrays
      for i,r in enumerate(res) :
        results[i][ri] = r[0]
        
    # Write out statistics over repeats
    print "\n%18s\t%8s\t%8s"%("","Tot-Mean","Stderr")
    for label,r in zip(args.labels,results) :
      print "%-18s\t%8.3f\t%8.3f"%(label,r.mean(),r.std()/np.sqrt(r.shape[0]))  

  
  