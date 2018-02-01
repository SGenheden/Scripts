# Author: Samuel Genheden samuel.genheden@gmail.com

import sys

import numpy as np

from sgenlib import parsing

if __name__ == '__main__':

    data = []
    for filename in sys.argv[1:] :
        data.append(parsing.parse2ndarray(filename))
    data = np.asarray(data)
    av = data.mean(axis=0)
    std = data.std(axis=0)/np.sqrt(data.shape[0])
    for irow in range(data.shape[1]) :
        print data[0, irow, 0],"\t",
        print "\t".join("%.3f\t%.3f"%(a,s)
            for a, s in zip(av[irow, 1:], std[irow, 1:]))
