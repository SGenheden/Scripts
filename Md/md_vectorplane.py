# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to calculate angle between the a user defined vector
and the membrane plane

Examples:
    md_vectorplane.py -f sim.dcd -s ref.pdb
"""

from sgenlib import moldyn
from sgenlib import mdactions

if __name__ == '__main__' :

    processor = moldyn.TrajectoryProcessor("Calculate the angle between a vector and the membrane plane")
    analysis = mdactions.VectorPlaneAngleAnalysis(processor)
    processor.setup(printargs=True)
    processor.process()
