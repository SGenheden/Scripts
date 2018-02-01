# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to calculate the orientation of a solute
Calculate the orientation for each solute, does not average

Examples:
    md_solorient.py -f sim.dcd -s ref.pdb --mask "name C1" "name C2"
    md_solorient.py -f sim.dcd -s ref.pdb --mask "name C1" "name C2" --axis 1.0 0.3 0.5
"""

from sgenlib import moldyn
from sgenlib import mdactions

if __name__ == '__main__' :

    processor = moldyn.TrajectoryProcessor("Calculte the solute orientation")
    analysis = mdactions.SoluteOrientation(processor)
    processor.setup(printargs=True)
    processor.process()
