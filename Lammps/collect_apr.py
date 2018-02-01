# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to collect APR results for a list of output

Examples:
collect_apr.py list
"""

import argparse
import glob
import logging
import os

import numpy as np

from sgenlib import parsing
from sgenlib import series
from sgenlib import umbrella

def _pull_win(str) :

    a, b = str.split("_p")
    return int(b.split(".")[0])

def _attach_win(str) :

    a, b = str.split("_a")
    return int(b.split(".")[0])


def _read_colvars_in(filename) :
    lines = []
    with open(filename, 'r') as f :
        lines = f.readlines()

    angfac = (180.0 / np.pi) * (180.0 / np.pi)
    for i, line in enumerate(lines) :
        if line.find("colvars d1") > -1 :
            kr = 0.5 * float(lines[i+1].strip().split()[1])
            r0 = float(lines[i+2].strip().split()[1])
        elif line.find("colvars a1") > -1 :
            kt = 0.5 * angfac * float(lines[i+1].strip().split()[1])
            t0 = float(lines[i+2].strip().split()[1])
        elif line.find("colvars a2") > -1 :
            kb = 0.5 * angfac * float(lines[i+1].strip().split()[1])
            b0 = float(lines[i+2].strip().split()[1])
    return kr, r0, kt, t0, kb, b0

def _read_datafile(filename, dt, skip) :
    """
    Read a data series and remove equilibration
    """
    data = parsing.parse2ndarray(filename)
    data[:,0] = data[:,0] * dt # Make the time column to ps
    freq = data[1,0] - data[0,0]
    nskip = int(skip/freq)
    d = data[(nskip+1):,:]
    return d

def _integrate_boot(windows, averages, errors, filename, logger, nboots=1000) :
    """
    Trapezoid integration with error estimation by bootstrapping
    """
    def _trapz(y, x) :
        return np.diff(x)*(y[:-1]+y[1:])*0.5

    logger.debug("Averages")
    for w, a, e in zip(windows, averages, errors) :
        logger.debug("%.3f %.3f %.3f"%(w,a,e))
    logger.debug("Maxium error %.3f for window %.3f"%(errors.max(),windows[np.argmax(errors)]))

    pmf = np.cumsum(_trapz(averages, windows))
    pmf_boots = np.zeros(nboots)
    err_boots = np.zeros((windows.shape[0]-1, nboots))
    for i in range(nboots) :
        av_boot = np.random.normal(averages, errors)
        err_boots[:,i] = _trapz(av_boot, windows)
        pmf_boots[i] = err_boots[:,i].sum()
    err_boots = np.sqrt(np.cumsum(err_boots.std(axis=1)**2))

    logger.debug("PMF")
    with open(filename, "w") as f :
        f.write("%.3f %.3f %.3f\n"%(windows[0],0,0))
        logger.debug("%.3f %.3f %.3f"%(windows[0],0,0))
        for w, p, e in zip(windows[1:], pmf, err_boots) :
            f.write("%.3f %.3f %.3f\n"%(w,p,e))
            logger.debug("%.3f %.3f %.3f"%(w,p,e))
    return pmf[-1], pmf_boots.std()

def _integrate_attach(filenames, system, dt, skip, logger):
    """
    Read an energy output from attach windows and integrate them
    """
    averages = [0]
    errors = [0]
    for filename in filenames :
        data = _read_datafile(filename, dt, skip)[:,1:].sum(axis=1)
        averages.append(data.mean())
        errors.append(series.standard_error(data))
    averages = np.asarray(averages)
    errors = np.asarray(errors)
    wins = [_attach_win(filename) for filename in filenames]
    wins.insert(0,0)
    wins = np.asarray(wins)*0.01
    logger.debug("\n--Attaching--")
    return _integrate_boot(wins, averages, errors, "%s_attach_pmf.dat"%system, logger)

def _integrate_pull(inputs, outputs, system, dt, skip, logger):
    """
    Read an energy output from pull windows and integrate them
    """
    logger.debug("\n--Pulling--")

    centers = []
    forces = []
    # Read input files to obtain equilibrium values and force constants
    for filename in inputs :
        kr, r0, kt, t0, kb, b0 = _read_colvars_in(filename)
        centers.append(r0)
        forces.append(kr)
    centers = np.asarray(centers)
    forces = np.asarray(forces)
    averages = []
    errors = []
    umbrellas = umbrella.UmbrellaSimulations(temperature=298)
    # Compute the work values from the output files
    for filename, center, force in zip(outputs, centers, forces) :
        data = -2*force*(_read_datafile(filename, dt, skip)[:,1]-center)
        averages.append(data.mean())
        errors.append(series.standard_error(data))
        umb = parsing.LammpsUmbrellaResults(center,weight=force*2,filename=filename)
        umb.skip(umb.samples.shape[0]-data.shape[0]-1,ispart=False)
        umbrellas.add(umb)
    averages = np.asarray(averages)
    errors = np.asarray(errors)

    umbrellas.make_histograms(nbins=100)
    umbrellas.plot_histograms(filename="%s_pull_hist.png"%system)
    overlaps = umbrellas.pairwise_overlap()
    logger.debug("Pair-wise overlap")
    for w1, w2, o in zip(centers[:-1],centers[1:],overlaps) :
        logger.debug("%.3f %.3f %.3f"%(w1, w2, o))
    mini = np.argmin(overlaps)
    logger.debug("Min overlap is %.3f for windows %d %d"%(overlaps.min(),centers[mini],

                                                            centers[mini+1]))
    z, pmf = umbrella.ExternalWham(umbrellas).pmf()
    ti, tierr =  _integrate_boot(centers, averages, errors, "%s_pull_pmf.dat"%system, logger)
    logger.debug("TI = %.3f WHAM = %.3f"%(ti,pmf[-1]-pmf[0]))
    return ti, tierr

def WRelease(kr,r0,kt,t0,kb,b0,T):
    """
       Do the analytical "release" of guest to standard concentration
    """
    ### SETUP
    # Example: print WRelease(5.0,24.00,100.0,180.0,100.0,180.0,298.15)
    # kr, r0: distance restraint (r) force constant, target value
    # kt, t0: angle restraint (theta) force constant, target value
    # kb, b0: angle restraint (b) force constant, target value
    # T: Temperature
    R = 1.987204118e-3 # kcal/mol-K, a.k.a. boltzman constant
    beta = 1/(T*R)
    rlb,rub,rst = [0.0,100.0,0.0001]  # r lower bound, upper bound, step size
    tlb,tub,tst = [0.0,np.pi,0.00005] # theta ",       ",            "
    blb,bub,bst = [0.0,np.pi,0.00005] # b     ",       ",            "
    def fr(val):
      return (val**2)*np.exp(-beta*kr*(val-r0)**2)
    def ft(val):
      return np.sin(val)*np.exp(-beta*kt*(val-np.radians(t0))**2)
    def fb(val):
      return np.sin(val)*np.exp(-beta*kb*(val-np.radians(b0))**2)
    ### Integrate
    rint,tint,bint = [0.0,0.0,0.0]
    intrange = np.arange(rlb,rub,rst)
    rint = np.trapz(fr(intrange),intrange)
    intrange = np.arange(tlb,tub,tst)
    tint = np.trapz(ft(intrange),intrange)
    intrange = np.arange(blb,bub,bst)
    bint = np.trapz(fb(intrange),intrange)
    return R*T*np.log(np.pi*(1.0/1660.0)*rint*tint*bint)


if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description="Script to collect output from APR")
    argparser.add_argument('outlist', help="the output list")
    argparser.add_argument('-outdir','--outdir',help="the directory with output files",default=".")
    argparser.add_argument('--dt', type=float, help="the timestep of the simulation", default=0.006)
    argparser.add_argument('--skip', type=float, help="the equilibration time", default=1200)
    argparser.add_argument('--temp', type=float, help="the temperature", default=298)
    args = argparser.parse_args()

    logger = logging.getLogger('collect_apr')
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(message)s')
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    console.setFormatter(formatter)
    logger.addHandler(console)
    logfile = logging.FileHandler("collect_apr.log", mode="w")
    logfile.setLevel(logging.DEBUG)
    logfile.setFormatter(formatter)
    logger.addHandler(logfile)

    systems = []
    with open(args.outlist, 'r') as f :
        systems = [line.strip() for line in f.readlines()]

    headers = "%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s"%("Attach","","Pull","","Release","Total","")
    print headers
    for sys in systems :
        logger.debug("#############\nSystem:\t%s"%sys)
        a_outputs = sorted(glob.glob(sys+"_rest_a*.out"), key=_attach_win)
        wattach, eattach = _integrate_attach(a_outputs, sys, args.dt, args.skip, logger)

        p_inputs = sorted(glob.glob(sys+"_colvars_p*.in"), key=_pull_win)
        p_outputs = sorted(glob.glob(sys+"_p*.colvars.traj"), key=_pull_win)
        wpull, epull = _integrate_pull(p_inputs, p_outputs, sys, args.dt, args.skip, logger)

        kr, r0, kt, t0, kb, b0 = _read_colvars_in(p_inputs[-1])
        wrelease = WRelease(kr, r0, kt, t0, kb, b0, args.temp)

        logger.debug("\n"+headers)
        logger.info("%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f"%(wattach,eattach,
                                        wpull,epull,wrelease,
                                        -(wattach+wpull+wrelease),np.sqrt(eattach**2+epull**2)))

        logger.debug("#############")
