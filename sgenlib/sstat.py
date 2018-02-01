# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Statistical routines
"""

import numpy as np
import scipy.stats as st

def repeated_oneway(data) :

  n = data.shape[0]
  k = data.shape[1]
  grand_mean = np.mean(data)
  measurement_mean = np.mean(data,axis=0)
  subject_mean = np.mean(data,axis=1)
  ssb = n*st.ss(measurement_mean-grand_mean)
#   ssw = st.ss(data-measurement_mean)
  ssw = np.sum(st.ss(data-measurement_mean))
  sss = k*st.ss(subject_mean-grand_mean)
  sse = ssw-sss
  dfb = k - 1
  dfe = (n-1)*(k-1)
  msb = ssb / float(dfb)
  mse = sse / float(dfe)
  f = msb / mse
  p = st.fprob(dfb,dfe,f)
  return f,p
  