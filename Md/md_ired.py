# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to calculate S2 order parameters using iRED

By default it calculates S2 of N-H vectors

Examples:
    md_ired.py -f sim.dcd -s ref.pdb
"""

from sgenlib import moldyn
from sgenlib import mdactions

if __name__ == '__main__' :

    processor = moldyn.TrajectoryProcessor("Calculate the iRED order parameters",
                    dosubsample=True)
    analysis = mdactions.IredAnalysis(processor)
    processor.setup()
    processor.process()
