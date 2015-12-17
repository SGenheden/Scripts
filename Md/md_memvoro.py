# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to perform Voronoi analysis of a membrane

At the moment it does two things:
  1) Calculate the area per lipid for each type of residue
  2) Count the number of neighbors for each pair of residue types

Examples
--------
md_memvoro.py -f md.xtc -s ref.gro --mask PO4 ROH --head PO4
"""

from sgenlib import moldyn
from sgenlib import mdactions

if __name__ == '__main__' :

    processor = moldyn.TrajectoryProcessor("Apply Voronoi analysis to membrane")
    analysis = mdactions.MemVoronoiAnalysis(processor)
    processor.setup(printargs=True)
    processor.process()
