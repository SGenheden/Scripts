# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to calculate the radius of gyration for a solute
Calculate the radius for each residue of the selection

Examples:
    md_gyration.py -f sim.dcd -s ref.pdb --mask "resname SOL"
"""

from sgenlib import moldyn
from sgenlib import mdactions

if __name__ == '__main__' :

    processor = moldyn.TrajectoryProcessor("Calculte the radius of gyration")
    analysis = mdactions.SoluteGyration(processor)
    processor.setup(printargs=True)
    processor.process()
