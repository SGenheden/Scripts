# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to calculate the RMSD of an atom selection

Examples:
md_rmsd.py -f md2_whole.xtc -s md1.gro
md_rmsd.py -f md2_whole.xtc -s md1.gro  --sel "protein and name CA"
"""

from sgenlib import moldyn
from sgenlib import mdactions

if __name__ == '__main__' :

    processor = moldyn.TrajectoryProcessor("Calculate RMSD of an atom selection")
    analysis = mdactions.RmsdAnalysis(processor)
    processor.setup(printargs=True)
    processor.process()
