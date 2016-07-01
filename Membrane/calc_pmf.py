# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to compute the potential of mean force (PMF) from
a set of umbrella simulations of small solutes in a membrane
"""

import argparse

import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
import numpy as np
import numpy.random as random
import scipy.stats as stats

import wham

if __name__ == '__main__' :

    parser = argparse.ArgumentParser(description="Calculationg PMF from umbrella sampling simulations")
    parser.add_argument('-f','--files',nargs="+",help="the output from the simulations.",default=[])
    parser.add_argument('-v0','--v0',nargs="+",help="the output from the simulations.",default=[])
    parser.add_argument('-v1','--v1',nargs="+",help="the output from the simulations.",default=[])
    parser.add_argument('-c','--centers',type=float,nargs="+",help="the centers of the simulations in Angstromgs",default=[])
    parser.add_argument('-m','--md',choices=["gmx","lammps","gmx_plumed"],help="the MD engine, should be either 'gmx' or 'lammps'",default="gmx")
    parser.add_argument('-o','--out',help="the prefix for the output",default="pmf")
    parser.add_argument('-w','--weight',type=float,help="the force constant used in the simulation",default=0.0)
    parser.add_argument('-t','--temp',type=float,help="the temperature used in the simulation",default=300.0)
    parser.add_argument('-r','--range',type=float,nargs="+",help="the range of the histograms")
    parser.add_argument('-n','--nbins',type=int,help="the number of bins of the histogram",default=200)
    parser.add_argument('-b','--blocks',type=int,help="the number of blocks used for averaging, default=5",default=5)
    parser.add_argument('-s','--skip',type=int,help="the number of snapshots to skip, default=0",default=0)
    parser.add_argument('--shorten',type=int,help="the number of snapshots to discard at the end, default=0",default=0)
    parser.add_argument('-d','--stride',type=int,help="the number of strides for the samples files, default=-1",default=-1)
    parser.add_argument('--synthetic',action="store_true",help="Use synthetic data rather than actual simulation data",default=False)
    parser.add_argument('--rehisto',action="store_true",help="Re-weight histograms rather than PMFs",default=False)
    parser.add_argument('--double',action="store_true",help="Indicates that there is two solutes in the membrane",default=False)
    parser.add_argument('--zconst',action="store_true",help="Turn on calculation with z-contraint method",default=False)
    parser.add_argument('--expand',action="store_true",help="Indicates that there is two solutes in the membrane and they should be concatenated",default=False)
    parser.add_argument('--diffusion',action="store_true",help="Turn on calculation of diffusion",default=False)
    parser.add_argument('--repeats',nargs="+",help="replacement pattern for multiple repeats")
    parser.add_argument('--dt',type=float,help="the time in ps per snapshot",default=10)
    args = parser.parse_args()

    if len(args.files) == 0 :
        print "No files specified so nothing to do."
        quit()

    if args.zconst :
        print "Sorry, but the z-contraint method is under developement!"
        quit()

    if args.zconst and not args.md =="gmx" :
        print "Z-constraint method only available with Gromacs outputfiles"
        quit()

    weight = args.weight
    if args.md in ["gmx_plumed","gmx"] : # Convert Gromacs weight to kcal/mol/A2
        weight = weight/4.184/100.0

    nskip = args.skip
    if args.skip == 0 and args.blocks > 1 :
        nskip = args.blocks+1

    # Setup constants
    boltzmann = 0.001982923700
    RT = boltzmann*args.temp

    # Setup stuff for double solutes
    if args.double :
        files2  = args.files[::-1]
        args.blocks = args.blocks*2

    # Setup simulation objects
    simulations = [wham.UmbrellaSimulations(args.temp) for b in range(args.blocks)]
    full_simulations = wham.UmbrellaSimulations(args.temp)

    # Decide which results file reader to use
    if args.md == "gmx" :
        results_class = wham.GromacsResultsFile
    elif args.md == "gmx_plumed" :
        results_class = wham.PlumedResultsFile
    elif args.md == "lammps" :
        results_class = wham.LammpsResultsFile

    if args.expand:
        for i, (filename, center) in enumerate(zip(args.files,args.centers[::-1])) :
            results = results_class(-center, weight, filename=filename, stride=args.stride, colidx=3, expansion=True, isforces=args.zconst)
            if args.shorten > 0 : results.shorten(args.shorten)
            results.skip(nskip)
            for sim, block in zip(simulations,results.block_it(args.blocks)) :
                sim.add(block)
            full_simulations.add(results)

    # Split distance files into blocks
    for i, (filename, center) in enumerate(zip(args.files,args.centers)) :
        results = results_class(center,weight,filename=filename,stride=args.stride,isforces=args.zconst)
        if args.shorten > 0 : results.shorten(args.shorten)
        results.skip(nskip)
        if args.repeats is not None :
            for rep in args.repeats[1:] :
                results2 = results_class(center,weight,filename=filename.replace(args.repeats[0],rep),stride=args.stride,isforces=args.zconst)
                if args.shorten > 0 : results2.shorten(args.shorten)
                results2.skip(nskip)
                results = results + results2

        if args.double :
            results2 = results_class(center,weight,filename=files2[i],stride=args.stride,colidx=3,isforces=args.zconst)
            if args.shorten > 0 : results2.shorten(args.shorten)
            results2.skip(nskip)
            if args.repeats is not None :
                for rep in args.repeats[1:] :
                    results3 = results_class(center,weight,filename=files2[i].replace(args.repeats[0],rep),stride=args.stride,colidx=3,isforces=args.zconst)
                    if args.shorten > 0 : results3.shorten(args.shorten)
                    results3.skip(nskip)
                    results2 = results2 + results3
            results = results + results2
        if args.synthetic : results.synthesize()

        if len(args.v0) == len(args.files) and len(args.v1) == len(args.files) :
            v0_data = wham.SimpleEnergyFile(center,weight,filename=args.v0[i])
            v1_data = wham.SimpleEnergyFile(center,weight,filename=args.v1[i])
            v0_data.samples = v1_data.samples-v0_data.samples
            v0_data.skip(nskip)
            print "Free energy difference = %.3f"%(-RT*np.log(np.mean(np.exp(-v0_data.samples/RT))))
            if args.rehisto :
                ene_blocks = v0_data.block_it(args.blocks)
            else :
                ene_blocks = [None for b in range(args.blocks)]
        else :
            ene_blocks = [None for b in range(args.blocks)]

        for sim,block,ene_block in zip(simulations,results.block_it(args.blocks),ene_blocks) :
            sim.add(block,energies=ene_block)

        full_simulations.add(results)
        if i == 0 :
            print "Will do calculations on %d snapshots / simulation"%(results.samples.shape[0])

    if not args.zconst :
        # Histogram the full data and calculate overlap
        full_simulations.make_histograms(args.nbins)
        full_simulations.plot_histograms(filename=args.out+"_hist.png")
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
            sim.make_histograms(args.nbins,boundaries=boundaries,weighted=args.rehisto)
            #sim.plot_histograms(filename=args.out+"_hist%d.png"%(i+1))
            if np.argmin(overlap) < len(args.files):
                j = np.argmin(overlap)
            else:
                j = np.argmin(overlap)-len(args.files)
            print "Minimum block overlap is %.2f between %s and %s"%(overlap.min(),args.files[j],args.files[j+1])
        pmf = wham.UmbrellaPmf()
        pmf.make(simulations,wham.ExternalWham)
        pmf.change_unit(wham.KJMOL)
        stride = 10

        if args.diffusion :
          diff = wham.UmbrellaPmf()
          diff.make(simulations,wham.Diffusion,no_offset=True)
          diff.av = diff.av  * 1E-4 / 1E-5
          diff.std = diff.std  * 1E-4 / 1E-5
          diff.write(args.out+"_diff.dat")
          diff.plot(filename=args.out+"_diff.png",stride=1,ylabel="D [10^-5 cm^2/s]")
    else :
        stride = 1
        pmf = wham.UmbrellaPmf()
        pmf.make(simulations,wham.ZConst,compute="pmf")
        pmf.change_unit(wham.KJMOL)

        resist = wham.UmbrellaPmf()
        resist.make(simulations,wham.ZConst,compute="resistance")
        resist.write(args.out+"_resist.dat")
        resist.plot(filename=args.out+"_resist.png",stride=stride)

    if args.md in ["gmx_plumed","gmx"] :
        filename = args.out+"_pmf.dat"
    elif args.md == "lammps" :
        filename = "wham."+args.out
    pmf.write(filename)
    pmf.plot(filename=args.out+"_pmf.png",stride=stride)
    print "Transfer free energy: %.3f +- %.3f"%pmf.transfer_dg()
    print "Water/lipid barrier: %.3f +- %.3f"%pmf.waterlipid_barrier()
    print "Penetration barrirer: %.3f +- %.3f"%pmf.penetration_barrier()
    print "Binding free energy: %.3f +- %.3f"%pmf.standard_dg()

    if args.double and args.blocks > 1 :
        half = args.blocks / 2
        pmf.average(end=half)
        print "First half:"
        print "\tTransfer free energy: %.3f +- %.3f"%pmf.transfer_dg()
        print "\tWater/lipid barrier: %.3f +- %.3f"%pmf.waterlipid_barrier()
        print "\tPenetration barrirer: %.3f +- %.3f"%pmf.penetration_barrier()
        print "\tBinding free energy: %.3f +- %.3f"%pmf.standard_dg()
        pmf.average(start=half)
        print "Second half:"
        print "\tTransfer free energy: %.3f +- %.3f"%pmf.transfer_dg()
        print "\tWater/lipid barrier: %.3f +- %.3f"%pmf.waterlipid_barrier()
        print "\tPenetration barrirer: %.3f +- %.3f"%pmf.penetration_barrier()
        print "\tBinding free energy: %.3f +- %.3f"%pmf.standard_dg()
