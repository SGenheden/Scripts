# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Routines to parse files and streams
"""

import fileinput
import re

import numpy as np


def stdin2ndarray() :
    """
    Parse columns and rows of data from standard input into a numpy.ndarray
    """
    data = []
    for line in fileinput.input() :
        if line[0] not in ["#","@"] :
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
    # Check for missing data at the end due to error
    cols0 = len(data[0])
    data = [line for line in data if len(line) == cols0]
    return np.array(data,dtype=float)

def parse_densityfile(filename) :
    """
    Parse densities from g_density
    """
    lines = []
    with open(filename, 'r') as f:
        lines = f.readlines()

    i = 0
    while lines[i].startswith("#") : i = i + 1

    legends = []
    while lines[i].startswith('@'):
        if lines[i].find("yaxis  label") > 0 :
            cols = lines[i].split('"')
            ylabel = cols[-2]
        elif re.search(" s\d+ legend",lines[i]):
            cols = lines[i].split('"')
            legends.append(cols[-2])
        i = i + 1

    data = []
    while i < len(lines):
        data.append(lines[i].strip().split())
        i = i + 1
    data = np.array(data, dtype=float)

    densities = {l : col for l, col in zip(legends, data[:,1:].T)}
    return data[:,0], densities, ylabel
