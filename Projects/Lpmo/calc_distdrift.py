# Author: Samuel Genheden samuel.genheden@gmail.com

import sys

import numpy as np

from sgenlib import parsing

data = parsing.parse2ndarray(sys.argv[1])
(nrows, nitems) = data.shape
nhalf = int(0.5*nrows)
for i in range(1, nitems):
    drift = np.abs(data[:nhalf, i].mean() - data[nhalf:, i].mean())
    print "%.3f\t%.3f"%(data[:, i].mean(), drift)
