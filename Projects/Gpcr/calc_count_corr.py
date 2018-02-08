# Author: Samuel Genheden samuel.genheden@gmail.com

import sys

import numpy as np
import scipy.stats as stats

from sgenlib import parsing

data = parsing.parse2ndarray(sys.argv[1])

r, prob = stats.pearsonr(60-data[:,0], data[:,1])
print "B2_i : %.2f %.2f"%(r, prob)
