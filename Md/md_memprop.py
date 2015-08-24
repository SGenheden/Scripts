# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to calculate APL, VPL and Dhh from an MD trajectory

Prints out the results in SI units

Examples:
    md_memprop.py -f md2.xtc -s md2.gro
    md_memprop.py -f md2.xtc -s md2.gro --watmask "resname W" --pmask "name PO4" --watvol 0.12
"""

from sgenlib import moldyn
from sgenlib import mdactions

if __name__ == '__main__' :

    processor = moldyn.TrajectoryProcessor("Calculate membrane properties")
    analysis = mdactions.MempropAnalysis(processor)
    processor.setup()
    processor.process()
