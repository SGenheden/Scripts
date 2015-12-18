# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Routines to perform analysis on series
"""

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
