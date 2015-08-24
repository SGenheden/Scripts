# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to calculate angle between the principal axis
of a selection and an arbitrary vector

By default is calculates the principal axis of the CA atoms
towards the z-axis.

Examples:
    md_principal.py -f sim.dcd -s ref.pdb
"""

from sgenlib import moldyn
from sgenlib import mdactions

if __name__ == '__main__' :

    processor = moldyn.TrajectoryProcessor("Calculte the principcal axis of a molecule and it angles with an axis")
    analysis = mdactions.PrincipalAxisAnalysis(processor)
    processor.setup(printargs=True)
    processor.process()
