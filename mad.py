# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program that prints the MAD of the first and second column from standard input.
The data on standard input should contain nothing but numbers in rows and columns.
"""
import sys

import numpy as np

from sgenlib import parsing

if __name__ == '__main__':

    domax = len(sys.argv) == 2 and sys.argv[1] == "x"
    domed = len(sys.argv) == 2 and sys.argv[1] == "m"
    sys.argv = [sys.argv[0]]
    data = parsing.stdin2ndarray()
    if not domax and not domed:
        print np.abs(data[:,0]-data[:,1]).mean()
    elif domax:
        print np.abs(data[:,0]-data[:,1]).max()
    elif domed:
        print np.median(np.abs(data[:,0]-data[:,1]))
