# Author: Samuel Genheden samuel.genheden@gmail.com

import argparse

import numpy as np

from sgenlib import umbrella
from sgenlib import parsing
from sgenlib.units import *

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Calculationg PMF from umbrella sampling simulations")
    parser.add_argument('-f','--files',nargs="+",help="the output from the simulations.",default=[])
    parser.add_argument('-c','--centers',type=float,nargs="+",help="the centers of the simulations in Angstromgs",default=[])
    parser.add_argument('-m','--md',choices=["gmx","lammps","gmx_plumed"],help="the MD engine, should be either 'gmx' or 'lammps'",default="gmx")
    parser.add_argument('-o','--out',help="the prefix for the output",default="pmf")
    parser.add_argument('-w','--weight',nargs="+", type=float,help="the force constant used in the simulation",default=[0.0])
    parser.add_argument('--temp',type=float,help="the temperature used in the simulation",default=300.0)
    parser.add_argument('--range',type=float,nargs="+",help="the range of the histograms")
    parser.add_argument('--nbins',type=int,help="the number of bins of the histogram",default=200)
    parser.add_argument('--blocks',type=int,help="the number of blocks used for averaging, default=5",default=5)
    parser.add_argument('--skip',type=int,help="the number of snapshots to skip, default=0",default=0)
    parser.add_argument('--shorten',type=int,help="the number of snapshots to discard at the end, default=0",default=0)
    parser.add_argument('--stride',type=int,help="the number of strides for the samples files, default=-1",default=-1)
    parser.add_argument('--repeats',nargs="+",help="replacement pattern for multiple repeats")
    args = parser.parse_args()

    if len(args.files) == 0 :
        print "No files specified so nothing to do."
        quit()


    if len(args.weight) == 1 :
        weights = np.asarray([args.weight]*len(args.files))
    else :
        weights = np.asarray(args.weight)
    if args.md in ["gmx_plumed","gmx"] : # Convert Gromacs weight to kcal/mol/A2
        weights = weights/4.184/100.0

    nskip = args.skip
    if args.skip == 0 and args.blocks > 1 :
        nskip = args.blocks+1

    # Setup simulation objects
    simulations = [umbrella.UmbrellaSimulations(args.temp) for b in range(args.blocks)]
    full_simulations = umbrella.UmbrellaSimulations(args.temp)

    # Decide which results file reader to use
    if args.md == "gmx" :
        results_class = parsing.GromacsUmbrellaResults
    elif args.md == "gmx_plumed" :
        results_class = parsing.PlumedUmbrellaResults
    elif args.md == "lammps" :
        results_class = parsing.LammpsUmbrellaResults


    # Split distance files into blocks
    for i, (filename, center, weight) in enumerate(zip(args.files,args.centers, weights)) :
        results = results_class(center,weight,filename=filename,stride=args.stride)
        if args.shorten > 0 : results.shorten(args.shorten)
        results.skip(nskip)
        if args.repeats is not None :
            for rep in args.repeats[1:] :
                results2 = results_class(center,weight,filename=filename.replace(args.repeats[0],rep),stride=args.stride)
                if args.shorten > 0 : results2.shorten(args.shorten)
                results2.skip(nskip)
                results = results + results2

        for sim, block in zip(simulations,results.block_it(args.blocks)) :
            sim.add(block)

        full_simulations.add(results)
        if i == 0 :
            print "Will do calculations on %d snapshots / simulation"%(results.samples.shape[0])

    # Histogram the full data and calculate overlap
    full_simulations.make_histograms(args.nbins)
    full_simulations.plot_histograms(filename=args.out+"_hist.png", xlabel="CV [A]")
    overlap = full_simulations.pairwise_overlap()
    if np.argmin(overlap) < len(args.files):
        j = np.argmin(overlap)
    else:
        j = np.argmin(overlap)-len(args.files)
    print "Minimum overlap is %.2f between %s and %s"%(overlap.min(),args.files[j],args.files[j+1])
    print "Minz = %.3f Maxz = %.3f Nbins = %d"%(full_simulations.bins[0],full_simulations.bins[-1],args.nbins)

    # Make histograms and run wham for each block
    if args.range is None :
        boundaries = [full_simulations.bins[0],full_simulations.bins[-1]]
    else :
        boundaries = args.range[:2]
    for i,sim in enumerate(simulations) :
        sim.make_histograms(args.nbins,boundaries=boundaries)
        if np.argmin(overlap) < len(args.files):
            j = np.argmin(overlap)
        else:
            j = np.argmin(overlap)-len(args.files)
        print "Minimum block overlap is %.2f between %s and %s"%(overlap.min(),args.files[j],args.files[j+1])
    pmf = umbrella.UmbrellaPmf()
    pmf.make(simulations,umbrella.ExternalWham)
    pmf.change_unit(KJMOL)
    stride = 10

    if args.md in ["gmx_plumed","gmx"] :
        filename = args.out+"_pmf.dat"
    elif args.md == "lammps" :
        filename = "wham."+args.out
    pmf.write(filename)
    pmf.plot(filename=args.out+"_pmf.png",stride=stride, xlabel="CV [A]")
