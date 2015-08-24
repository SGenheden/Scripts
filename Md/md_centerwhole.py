# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to transform a molecular dynamics trajectory.

It can
1) Make molecules whole over periodic boxes
2) Center a protein in the middle of the box
3) Align protein with an RMSD fit

Works only for rectangular boxes.

Examples:
    md_centerwhole.py -f sim.dcd -s ref.pdb
    md_centerwhole.py -f sim.dcd -s ref.pdb --noalign
"""

from sgenlib import moldyn
from sgenlib import mdactions

if __name__ == '__main__' :

  processor = moldyn.TrajectoryProcessor("Center and make a trajectory whole")
  analysis = mdactions.CenterWholeAlign(processor)
  processor.setup(printargs=True)
  processor.process()
