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
import sys
import argparse

import numpy as np
import matplotlib
if not "DISPLAY" in os.environ or os.environ["DISPLAY"] == "" :
  matplotlib.use('Agg')
import matplotlib.pyplot as plt

import gpcr_lib
from sgenlib import colors

def moving(data,window) :
  weights = np.repeat(1.0,window) / float(window)
  return np.convolve(data,weights,'valid')

def _anal_contacts(filenames,labels,plot=[],time=None,window=None) :
  """
  Main analysis routine
  
  Parameters
  ----------
  filenames : list of strings
    the files to read and analyse
  labels : list of strings
    the labels for the files
  plot : list of integers, optional
    the series to plot
  time : float, optional
    the total simulation time
  window : float, optional
    the window time for averaging
    
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
  for i,(state,label) in enumerate(zip(states,labels),1) : 
    counton = state.sum(axis=1)
    print "%-18s\t%8.3f\t%8.3f"%(label,counton.mean(),counton.std())  
    results.append(np.array([counton.mean(),counton.std()]))
    
    if len(plot) > 0 and time is not None and i in plot :
      ns_per_snapshot = time / float(counton.shape[0])
      snapshot_per_window = window / ns_per_snapshot
      series = moving(counton,snapshot_per_window)
      f = plt.figure(i)
      x = np.arange(series.shape[0])*ns_per_snapshot
      f.gca().plot(x,series,color=colors.color(plot.index(i)))
      f.gca().set_xlabel("Time (ns)")
      f.gca().set_ylabel(label)
      f.savefig("series%d.png"%(i),format="png",dpi=300)
    
  return results

if __name__ == '__main__' :

  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Analysing counts from state files")
  parser.add_argument('-f','--files',nargs='+',help="a list of input files.",default=[])
  parser.add_argument('-l','--labels',nargs='+',help="a label for each state.",default=[])
  parser.add_argument('--repeats',nargs="+",help="replacement pattern for multiple repeats",default=["r1_","r2_","r3_","r4_","r5_"]) 
  parser.add_argument('--time',type=float,help="total simulation time in ns",default=50000) 
  parser.add_argument('--window',type=float,help="the window time in ns",default=50) 
  parser.add_argument('--plot',type=int,nargs="+",help="which series to plot",default=[])
  args = parser.parse_args()
  
  res0 = _anal_contacts(args.files,args.labels,args.plot,args.time,args.window)
  
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

  
  
