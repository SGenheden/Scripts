# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Routines to perform analysis on series
"""

import numpy as np
from scipy.stats import spearmanr

def find_equilibration(x,y,atleast=50,threshold=0.05,nperm=0) :
  """
  Find the equilibration time of a data series

  Parameters
  ----------
  x : Numpy array
    the x data
  y : Numpy array
    the y data
  atleast : int, optional
    this is the smallest number of snapshots considered to be production
  threshold : float, optional
    the confidence level
  nperm : int, optional
    if greater than zero, a permutation test is peformed with nperm synthetic permuations

  Returns
  -------
  int
    the length of the equilibration period
  """
  for i in range(0,y.shape[0]-atleast) :
    if nperm == 0 : # Performs an assymptotic test, fast
      tau,p = spearmanr(x[i:],y[i:])
      if p > threshold : return i
    else : # Performs a rigorous permutation test, slow
      rhos = np.zeros(nperm)
      rho0,p = spearmanr(x[i:],y[i:])
      for j in range(nperm) :
        rhos[j],p = spearmanr(x[i:],np.random.permutation(y[i:]))
      ncnt = (rho0*rho0 < rhos*rhos).sum()
      if float(ncnt) / float(nperm) > threshold : return i
  return y.shape[0]-atleast

def standard_error(y) :
    """
    Find the standard error of a time series by block averaging

    Parameters
    ----------
    y : ndarray
        the y data

    Returns
    -------
    float
      the standard error
    """
    def _factors(n):
        # Return list of integer factors
        factor = []
        sqrt = int(round(np.sqrt(n) + 0.5))
        i = 1
        while i <= sqrt:
            if n % i == 0:
                factor.append(i)
                j = n / i
                if j != i:
                    factor.append(j)
            i += 1
        return sorted(factor, key=int)

    n = y.shape[0]
    std_max = 0.0
    for fac in _factors(n)[:-2] :
        blocks = y.reshape((-1, fac))
        block_av =  blocks.mean(axis=1)
        std = np.std(block_av, ddof=0) / np.sqrt(blocks.shape[0] - 1)
        std_max = max(std, std_max)
    return std_max
