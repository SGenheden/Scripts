# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program that prints block error of each column of
data parsed from standard input. The data on standard input should contain
nothing but numbers in rows and columns. The number of block is two by
default but can be supplied at the command line
"""

import sys, math
from sgenlib import parsing

if __name__ == '__main__':

    stde = len(sys.argv) > 1 and sys.argv[1] == "-e"
    if stde :
        sys.argv.pop(1)
    nblocks = 2 if len(sys.argv) == 1 else int(sys.argv[1])
    sys.argv = [sys.argv[0]]
    data = parsing.stdin2ndarray()
    stride = data.shape[0] / nblocks
    data = data[stride::stride,:].T
    if stde :
        print "\t".join("%10.4f"%(col.std()/math.sqrt(len(col))) for col in data)
    else :
        print "\t".join("%10.4f"%(col.std()) for col in data)
