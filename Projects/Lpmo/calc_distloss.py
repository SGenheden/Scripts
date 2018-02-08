# Author: Samuel Genheden samuel.genheden@gmail.com

import sys

import numpy as np

from sgenlib import parsing

data = parsing.parse2ndarray(sys.argv[1])[:, 1:]
(nrows, nitems) = data.shape
nhalf = int(0.5*nrows)

diff = np.median(np.abs(data[0, : ] - data[:, :]), axis=1)
print diff.shape
for j, d in enumerate(diff) :
    if d > 0.1 :
        print "%0.3f %d"%(1+j*0.01, j)
        print data[j, :]
        break
