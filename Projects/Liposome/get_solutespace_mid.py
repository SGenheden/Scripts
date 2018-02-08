# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Script to find the space where to insert the solutes in the middle of the liposome
"""

import sys

import MDAnalysis
import MDAnalysis.lib.distances as mddist
import numpy as np

u = MDAnalysis.Universe(sys.argv[1])
inner_lipids = u.select_atoms("name PO4 and resid 8161:9107" )
outer_lipids = u.select_atoms("name PO4 and resid 9108:10688")
inner_com = np.asarray([inner_lipids.center_of_geometry()])
outer_com = np.asarray([outer_lipids.center_of_geometry()])
inner_radius = mddist.distance_array(inner_com,inner_lipids.positions,None).mean()
outer_radius = mddist.distance_array(outer_com,outer_lipids.positions,None).mean()
mid = inner_radius + 0.5*(outer_radius - inner_radius)

print "outside sphere %.3f %.3f %.3f %.3f"%(outer_com[0,0], outer_com[0,1], outer_com[0,2], mid-2.5)
print "inside sphere %.3f %.3f %.3f %.3f"%(outer_com[0,0], outer_com[0,1], outer_com[0,2], mid+2.5)
