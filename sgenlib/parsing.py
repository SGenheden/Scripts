# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Routines to parse files and streams
"""

import fileinput

import numpy as np


def stdin2ndarray() :
    """
    Parse columns and rows of data from standard input into a numpy.ndarray
    """
    data = []
    for line in fileinput.input() :
        data.append(line.strip().split())

    return np.array(data,dtype=float)

def parse2ndarray(filename):
    """
    Parse columns and rows of data from a file into a numpy.ndarray

    Treats lines starting with # and @ as comments
    """
    data = []
    with open(filename,"r") as f :
        data = [line.strip().split() for line in f.readlines()
                 if line[0] not in ["#","@"] ]
    return np.array(data,dtype=float)
