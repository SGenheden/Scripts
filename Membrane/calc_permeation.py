import argparse

import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
import numpy as np
import scipy.stats as stats
import numpy.random as random

import wham

if __name__ == '__main__' :

  # Command-line input

  parser = argparse.ArgumentParser(description="Calculationg permeation from PMF and diffusion profile")
  parser.add_argument('-p','--pmf',help="the PMF from the simulations.")
  parser.add_argument('-d','--diff',help="the diffusion from the simulations.")
  parser.add_argument('-o','--out',help="the output prefix",default="")
  args = parser.parse_args()

  pmf = wham.UmbrellaPmf()
  pmf.read(args.pmf)

  diff = wham.UmbrellaPmf()
  diff.read(args.diff)
  diff.av = diff.av*1E-5
  diff.std = diff.std*1E-5

  fpmf = plt.figure(1)
  fpmf.gca().plot(pmf.z,pmf.av)
  fpmf.savefig(args.out+"perm_pmf.png")

  fdiff = plt.figure(2) 
  fdiff.gca().plot(diff.z,diff.av/1E-5)
  sel = np.logical_and(pmf.z>diff.z[0],pmf.z<diff.z[-1])
  from scipy.interpolate import interp1d
  f = interp1d(diff.z,diff.av,kind="cubic")
  diff2 = f(pmf.z[sel])
  fdiff.gca().plot(pmf.z[sel],diff2/1E-5)
  fdiff.savefig(args.out+"perm_diff.png")

  resistance = np.exp(pmf.av[sel]/pmf.kT) / diff2
  frest = plt.figure(3)
  frest.gca().plot(pmf.z[sel],resistance/1E6)
  frest.savefig(args.out+"perm_rest.png")

  print 1/np.trapz(resistance,pmf.z[sel])
