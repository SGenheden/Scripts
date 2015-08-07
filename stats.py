# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program that prints the average and standard deviation of each column of 
data parsed from standard input. The data on standard input should contain 
nothing but numbers in rows and columns.
"""

from sgenlib import parsing

if __name__ == '__main__':

    data = parsing.stdin2ndarray()
    print "\t".join("%10.4f %10.4f"%(col.mean(),col.std()) for col in data.T)
