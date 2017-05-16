# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to strip atoms from a trajectory

Examples:
md_strip.py -f md2_whole.xtc -s md1.gro --sel "name SOL"
"""

from sgenlib import moldyn
from sgenlib import mdactions

if __name__ == '__main__' :

    processor = moldyn.TrajectoryProcessor("Strip atoms from a trajectory")
    analysis = mdactions.StripAtoms(processor)
    processor.setup(printargs=True)
    processor.process()
