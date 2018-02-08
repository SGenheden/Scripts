# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to calculate the vdW energy between a ROH bead and protein beads
on a rectangular grid

Writes out the 3D grid in dx-format and an averaged 2D plot in png-format

Example:
  gpcr_plot_gridenergy.py -p system.top -m b2
"""

import argparse
import os
import sys

import numpy as np
import matplotlib
if not "DISPLAY" in os.environ or os.environ["DISPLAY"] == "" :
  matplotlib.use('Agg')
import matplotlib.pylab as plt

import  gpcr_lib
from sgenlib import gmx

def _writeDX(grid,origin,spacing,filename) :
  f = open(filename, 'w')
  f.write("object 1 class gridpositions counts %5d%5d%5d\n"%(grid.shape[0],grid.shape[1],grid.shape[2]))
  f.write("origin %9.4f%9.4f%9.4f\n"%(origin[0],origin[1],origin[2]))
  f.write("delta %10.7f 0.0 0.0\n"%spacing)
  f.write("delta 0.0 %10.7f 0.0\n"%spacing)
  f.write("delta 0.0 0.0 %10.7f\n"%spacing)
  f.write("object 2 class gridconnections counts %5d%5d%5d\n"%(grid.shape[0],grid.shape[1],grid.shape[2]))
  f.write("object 3 class array type double rank 0 items  %10d data follows\n"%(grid.shape[0]*grid.shape[1]*grid.shape[2]))
  cnt = 0
  for x in range(grid.shape[0]) :
    for y in range(grid.shape[1]) :
      for z in range(grid.shape[2]) :
        f.write("%19.10E"%grid[x,y,z])
        cnt = cnt + 1
        if cnt >= 3 :
          cnt = 0
          f.write("\n")
  if cnt > 0 : f.write("\n")
  f.write('attribute "dep" string "positions"\n')
  f.write('object "regular positions regular connections" class field\n')
  f.write('component "positions" value 1\n')
  f.write('component "connections" value 2\n')
  f.write('component "data" value 3\n')
  f.close()

if __name__ == '__main__' :

  # Command-line input
  parser = argparse.ArgumentParser(description="Analyzing GPCR grid interactions")
  parser.add_argument('-p','--top',help="a topology file",default="system.top")
  parser.add_argument('-m','--mol',choices=["b2","a2a","b2_a","a2a_a"],help="the protein molecule")
  args = parser.parse_args()

  # Read a Xrya structure file
  xray = gpcr_lib.load_xray(args.mol,loadsigma=True)
  pdb = xray.pdbfile

  # Remove cholesterols at the end
  while pdb.residues[-1].resname == "CHOL" :
    for atom in pdb.residues[-1].atoms :
      pdb.atoms.remove(atom)
    del pdb.residues[-1]
  pdb.xyz = np.array([atom.xyz for atom in pdb.atoms])

  # Read a Gromacs topology file
  top = gmx.TopFile(args.top)

  # Store away all interesting LJ parameter pairs
  roh_type = "SP1"  #"SC1"
  roh_pairs = []
  for p in top.pairtypes :
    if p.atomi == roh_type or p.atomj == roh_type : roh_pairs.append(p)

  # Find the protein molecule definition
  protein_type = None
  for m in top.moleculetypes :
    if m.name == "Protein" : protein_type = m
  if protein_type is None :
    raise Exception("Could not find a protein molecule type")

  # Check that the given structure is compatible with the found protein force field
  if len(protein_type.atoms) != len(pdb.atoms) :
    raise Exception("Different number of atoms in ff (%d) and structure (%d)"%(len(protein_type.atoms),len(pdb.atoms)))

  # Store away LJ c6 and c12 coefficients
  c6 = np.zeros(len(protein_type.atoms))
  c12 = np.zeros(len(protein_type.atoms))
  sigmas = np.zeros(len(protein_type.atoms))
  for i,atom in enumerate(protein_type.atoms) :
    found = False
    for p in top.pairtypes :
      if p.atomi == atom.type and p.atomj == atom.type :
        sigmas[i] = (p.epsilon / p.sigma)**(1.0/6.0)*10.0/2.0
        break
    for p in roh_pairs :
      if p.atomi == atom.type or p.atomj == atom.type :
        found = True
        c6[i] = p.sigma
        c12[i] = p.epsilon
    if not found :
      print "Incomplete force field spec for type %s"%atom.type
      quit()

  # Store some coordinate stats
  mincoord = pdb.xyz.min(axis=0)
  maxcoord = pdb.xyz.max(axis=0)
  midcoord = pdb.xyz.mean(axis=0)

  # Determine the size of the box
  zlen = int(np.round(maxcoord[2]-mincoord[2]))
  if zlen % 2 != 0 : zlen = zlen + 1
  xylen = 70
  grid = np.zeros([xylen,xylen,zlen])

  # Convert coordinates to nm
  pdbxyz = pdb.xyz / 10.0

  # Loop over all the entire grid
  for x in range(-xylen/2,xylen/2+1) :
    for y in range(-xylen/2,xylen/2+1) :
      for z in range(-zlen/2,zlen/2+1) :
         gridcoord = (midcoord + np.array([x,y,z]))/10.0
         diff2 = np.sum((pdbxyz-gridcoord)**2,axis=1)
         diff6 = diff2*diff2*diff2
         diff12 = diff6*diff6
         lj = c12/diff12-c6/diff6
         lj[diff2>(1.2*1.2)] = 0.0
         grid[x+xylen/2-1,y+xylen/2-1,z+zlen/2-1] = lj.sum()
         #print gridcoord*10.0,lj.sum()

  grid[grid>0.0] = 0.0

  # Write out the 3D grid
  origin = midcoord + np.array([-xylen/2,-xylen/2,-zlen/2])
  _writeDX(grid,origin,1.0,"gridenergy.dx")

  # Plot 2D plot as average over z-dim
  grid_low = grid[:,:,:zlen/2].mean(axis=2)
  grid_upp = grid[:,:,zlen/2:].mean(axis=2)

  fig = plt.figure(1)
  for i,(grid,leaflet) in enumerate(zip([grid_low,grid_upp],["low","upp"]),1) :
    a = fig.add_subplot(1,2,i)
    gridrm = np.ones([70,70,4])
    gridrm[grid<0.0,3] = 0.0
    if leaflet == "upp" :
      grid = grid[:,::-1]
      gridrm = gridrm[:,::-1]
    im = a.imshow(grid,extent=[-35,35,-35,35],origin="lower",cmap=plt.cm.YlOrRd_r)
    a.imshow(gridrm,extent=[-35,35,-35,35],origin="lower")
    gpcr_lib.plot_density_xray(a,None,"",0,0,xray,leaflet,gpcr_lib.side_name[leaflet].capitalize(),number=None,plotn=False,drawchol=False)

  ax = gpcr_lib.draw_colormap(fig,im, text='')
  ax.text(1.03,0.80,'$\mathrm{Energy}$')

  fig.savefig("gridenergy.png",format="png")
