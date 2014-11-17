# Author: Samuel Genheden, 2012-2014

"""
Program to bootstrap quality metrics
"""

import argparse

import numpy as np

import quality

def _parse_input(infile,errfile,expectfile,uncert,verbose) :
  """
  Parse the user input into numpy arrays
  """
  p = []
  pe = []
  e = []
  ee = []

  try :
    lines = open(infile,"r").readlines()
  except :
    raise Exception("Error: Problem reading a file: %s"%infile)
  ncols = -1
  for s in lines :
    if s[0] != "#" :
      cols = s.strip().split()
      if ncols == -1 and len(cols) > 0:
        ncols = len(cols)
        if errfile == None and expectfile == None and ncols < 3 :
          raise Exception("Error: Too few columns in the input file. \nIf the -e and -x options are not set at the command-line, each line in the input file must contain at least 3 columns.")
        elif errfile == None and expectfile == None and ncols > 4 and verbose:
          print "Warning: Since the -e and -x options are not set, columns 5 and above will be ignored file %s."%infile
      if len(cols) > 0 and len(cols) != ncols and verbose:
        print "Warning: Skipping uncommented line of none-zero length (%s) in file %s"%(s,infile)
        print "         This line has not an qual number of space-separated columns (%d) than the first (%d)"%(len(cols),ncols)
      elif len(cols) > 0 :
        if errfile == None and expectfile == None :
          p.append(cols[0])
          pe.append(cols[1])
          e.append(cols[2])
          if ncols > 3:
            ee.append(cols[3])  
        elif errfile == None and expectfile != None :
          if (ncols % 2) != 0 :
            raise Exception("Error: Invalid number of columns in the input file. \nIf the -x option is set but not the -e option, each line in the input must have an even number of columns.")
          if ncols > 2 :
            p.append(list(cols[0:ncols:2]))
            pe.append(list(cols[1:ncols:2]))
          else :
            p.append(cols[0])
            pe.append(cols[1])
        elif errfile != None and expectfile != None :
          if ncols > 1 :
            p.append(list(cols))
          else :
            p.append(cols[0])

  if errfile != None :     
    pe = []
    try :
      lines = open(errfile,"r").readlines()
    except :
      raise Exception("Error: Problem reading a file: %s"%errfile)
    ncols = -1
    for s in lines :
      if s[0] != "#" :
        cols = s.strip().split()
        if ncols == -1 and len(cols) > 0 :
          ncols = len(cols)
          if ncols != len(p[0]) :
            raise Exception("Error: Invalid columns in file %s. \nThe number of columns (%d) must be equal to the number of columns in the input file (%d)"%(errfile,ncols,len(p[0])))
        if len(cols) > 0 and len(cols) != ncols and verbose :
          print "Warning: Skipping uncommented line of none-zero length (%s) in file %s"%(s,errfile)
          print "         This line has not an equal number of space-separated columns (%d) than the first (%d)"%(len(cols),ncols)
        elif len(cols) > 0 :
          pe.append(list(cols))
    if len(pe) != len(p) :
      raise Exception("Error: Invalid number of rows (%d) in %s. \nYou must supply error values for the exact number of items in the inputfile (%d)."%(len(pe),errfile,len(p)))
 
  if expectfile != None:
    e = []
    ee = []
    try :
      lines = open(expectfile,"r").readlines()
    except :
      raise Exception("Error: Problem reading a file: %s"%expectfile)
    ncols = -1
    for s in lines :
      if s[0] != "#" :
        cols = s.strip().split()
        if ncols == -1 and len(cols) > 0 :
          ncols = len(cols)
          if ncols > 2 and verbose :
            print "Warning: Columns 3 and above will be ignored file %s."%expectfile
        if len(cols) > 0 and len(cols) != ncols and verbose:
          print "Warning: Skipping uncommented line of none-zero length (%s) in file %s"%(s,expectfile)
          print "         This line has not an equal number of space-separated columns (%d) than the first (%d)"%(len(cols),ncols)
        elif len(cols) > 0 :
          e.append(cols[0])
          if len(cols) > 1 :
            ee.append(cols[0])
    if len(e) != len(p) :
      raise Exception("Error: Invalid number of rows (%d) in %s. \nYou must supply expected values for the exact number of items in the inputfile (%d)."%(len(e),expectfile,len(p)))
 
  if uncert != None :
    ee = np.ones(len(e))*uncert

  pred = np.array(p,dtype=float)
  prederr = np.array(pe,dtype=float)
  expect = np.array(e,dtype=float)
  expecterr = np.array(ee,dtype=float)

  try :
    ncols = pred.shape[1]
  except :
    ncols = 1

  if verbose :
    print "Read %d columns and %d rows of predicted values from %s"%(ncols,len(p),infile)
    if errfile == None :
      print "Read prediction errors from %s"%infile
    else :
      print "Read prediction errors from %s"%errfile
    if expectfile == None :
      print "Read expected values from %s"%infile
      if uncert == None and len(ee) > 0 :
        print "Read errors of expected values from %s"%infile
    else :
      print "Read expected values from %s"%expectfile
      if uncert == None and len(ee) > 0 :
        print "Read errors of expected values from %s"%expectfile
    if uncert != None :
      print "Set the errors of expected values to %.3f"%uncert
  
  return (pred,prederr,expect,expecterr,ncols)
  
if __name__ == '__main__':

  disclaimer="""
If -i and -m is not supplied on the command-line, the program will ask for them
\n
The program can take input from many sources, it can parse a single or a series 
of predictions.

 For single prediction
   - Predictions are specified with the -i option and parsed from the 1:st column
   - Prediction Error can either be specified with the -e option or parsed from
     the 2:nd column of the -i option
   - Expected values can either be specified with the -x option or parsed from 
     the 3:rd column of the -i option

 For series
  - Predictions are specified with -i option, and depending on if the -e option
    is used the entire file or 1:st, 3:rd, etc columns are parsed
  - Prediction Error can be specified with the -e option or with the -i option,
    in which case the errors are parsed from the 2:nd, 4:th, etc columns
  - Expected values must be specified with the -x option

Errors of expected values can either be specified as the 2:nd column in the -x
 option, or the 4:th column in the -i option (single prediction only). It can
 also be specified as a single value with the -u option, and this value is 
 then applied to all data points.

\n
If a question mark (?) is specified as a metric, all available metrics
 will be printed to screen.
  """
  # Parse the command-line
  parser = argparse.ArgumentParser(description="Program to calculate quality metrics",epilog=disclaimer,formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument('-i', '--input',help="the input file")
  parser.add_argument('-e', '--errors',help="a file with prediction errors")
  parser.add_argument('-x', '--expected',help="a file with expected values")
  parser.add_argument('-m', '--metric',nargs='+',help="a list of metrics to compute")
  parser.add_argument('-b', '--boots',type=int,help="number of bootstraps",default=1000)
  parser.add_argument('-s', '--scaling',type=int,help="scaling factor for prediction errors",default=1)
  parser.add_argument('-l', '--level',type=float,help="level of confidence",default=1.64500)
  parser.add_argument('-u', '--uncertainty',type=float,help="experimental uncertainty")
  parser.add_argument('-v', dest='verbose', action='store_true',help="turns on verbose output")
  parser.add_argument('-pd', dest='plotDelta', action='store_true',help="plot delta distribution (only for a single series)")
  parser.add_argument('-ps', dest='plotSample', action='store_true',help="plot sample distribution (only for a single series)")
  args = parser.parse_args()
  
  # Ask for the essential input options if they were not supplied on the command-line
  if args.input == None :
    print "Enter the input filename: ",
    args.input = raw_input()
  if args.metric == None :
    print "Enter the quality metrics to compute: ",
    args.metric = raw_input().strip().split()
  
  # Initialize the global dictionary
  metrics = []
  for m in args.metric :
    if m == "?" :
      print "\nAvailable metric keys are: ",quality.QualityCollection.metrickeys()
      quit()
    else :
      metrics.append(m.lower())
  
  # Parse input data to obtain predictions and expected values as well as their uncertainties
  try :
    (pred,prederr,expect,expecterr,nseries) = _parse_input(args.input,args.errors,args.expected,args.uncertainty,args.verbose)
  except Exception as e :
    print e.args[0]
    quit()

  # Optionally calculate standard errors from standard deviations
  if args.scaling > 1 :
    if args.verbose : print "Scaling all prediction errors with %.3f"%np.sqrt(args.scaling)
    prederr = prederr / np.sqrt(args.scaling)

  # Create QualityCollection object(s), do bootstrapping and print out the results
  if nseries == 1 :
    qc=quality.QualityCollection(pred,prederr,expect,expecterr,metrics,nboots=args.boots,verbose=args.verbose,level=args.level)
    qc.bootstrap()
    qc.printH0()
    qc.printStat()
    qc.plotDistributions(plotDelta=args.plotDelta,plotSample=args.plotSample)
  else :
    qc0=quality.QualityCollection(pred[:,0],prederr[:,0],expect,expecterr,metrics,nboots=args.boots,verbose=args.verbose,level=args.level)
    qc0.printH0()
    biased = np.zeros((len(qc0.metriclst),nseries))
    std = np.zeros(biased.shape)
    for s in range(nseries) :
      qc=quality.QualityCollection(pred[:,s],prederr[:,s],expect,expecterr,metrics,nboots=args.boots,verbose=args.verbose,level=args.level)
      qc.bootstrap()
       # Collect results for each of the metrics
      for m in range(len(qc.metriclst)) :
        biased[m,s] = qc.metriclst[m].biased
        std[m,s] = qc.metriclst[m].std
    # Print results
    print "%10s %9s"%("","Biased"),
    for s in range(nseries-1) :
      print " %9s"%(""),
    print " %9s"%"Stdev",
    for s in range(nseries-1) :
      print " %9s"%(""),
    print ""
    for m in range(biased.shape[0]) :
      print "%-10s"%qc0.metriclst[m].name(),
      for s in range(nseries) :
        print " %9.4f"%biased[m,s],
      for s in range(nseries) :
        print " %9.4f"%std[m,s],
      print ""
      