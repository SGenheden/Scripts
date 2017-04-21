# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to calculate membrane densities and properties from them

Prints out the results in SI units

Examples:
    md_memdens.py -f md2.xtc -s md2.gro
"""

from sgenlib import moldyn
from sgenlib import mdactions

if __name__ == '__main__' :

    processor = moldyn.TrajectoryProcessor("Calculate membrane densitiy properties",
        dosubsample=True)
    analysis = mdactions.MemDensAnalysis(processor)
    processor.setup()
    processor.process()
