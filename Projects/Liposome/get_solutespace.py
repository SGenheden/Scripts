# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Script to find the space where to insert the solutes
"""

import sys

import MDAnalysis
import MDAnalysis.lib.distances as mddist
import numpy as np

u = MDAnalysis.Universe(sys.argv[1])
lipids = u.select_atoms("name PO4 and resid 9108:10688")
com = np.asarray([lipids.center_of_geometry()])
radius = mddist.distance_array(com,lipids.positions,None).mean()

print "outside sphere %.3f %.3f %.3f %.3f"%(com[0,0], com[0,1], com[0,2], radius+10)
print "inside box 0.0 0.0 0.0 %.3f %.3f %.3f"%tuple(u.dimensions[:3])
