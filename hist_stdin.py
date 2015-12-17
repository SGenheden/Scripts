# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program that histograms the data parsed from standard input.
The data on standard input should contain nothing but numbers in a single column.
"""

import matplotlib.pyplot as plt

from sgenlib import parsing

if __name__ == '__main__':

    data = parsing.stdin2ndarray()
    nbins = int(float(data.shape[0])/100.0)
    plt.hist(data,bins=nbins)
    plt.show()
