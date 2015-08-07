# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program that prints the sum of each column of data parsed from standard input.
The data on standard input should contain nothing but numbers in rows and columns.
"""

from sgenlib import parsing

if __name__ == '__main__':

    data = parsing.stdin2ndarray()
    print "\t".join("%10.4f"%col.sum() for col in data.T)
