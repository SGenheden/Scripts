# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to calculate C-C chain order parameter from MD trajectory

The carbon atoms that form the tail on which the order parameters
should be calculated are given as command-line arguments. One or
more chains can be given.

Examples
--------
md_chainorder.py -f md2_whole.xtc -s md1.gro -c POPC:C1A-D2A-C3A-C4A POPC:C1B-C2B-C3B-C4B
"""

from sgenlib import moldyn
from sgenlib import mdactions

if __name__ == '__main__':

    processor = moldyn.TrajectoryProcessor("Calculate tail order parameter from MD trajectory")
    analysis = mdactions.ChainOrderAnalysis(processor)
    processor.setup(printargs=True)
    processor.process()
