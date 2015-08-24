# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to calculate RMSF of an MD trajectory

By defaults calculates the RMSF of the C, CA and N atoms

Examples:
    md_rmsf.py -f sim.dcd -s ref.pdb
"""

from sgenlib import moldyn
from sgenlib import mdactions

if __name__ == '__main__' :

    processor = moldyn.TrajectoryProcessor("Calculate the bfactor of atoms")
    analysis = mdactions.RMSFAnalysis(processor)
    processor.setup(printargs=True)
    processor.process()
