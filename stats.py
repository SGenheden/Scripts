# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program that prints the average and standard deviation of each column of
data parsed from standard input. The data on standard input should contain
nothing but numbers in rows and columns.
"""

import sys
from sgenlib import parsing

if __name__ == '__main__':

    invert = len(sys.argv) == 2 and sys.argv[1] == "r"
    sys.argv = [sys.argv[0]]
    data = parsing.stdin2ndarray()
    if not invert:
        data = data.T
        sep = "\t"
    else:
        sep = "\n"
    print sep.join("%10.4f\t%10.4f"%(col.mean(),col.std()) for col in data)
