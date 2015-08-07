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