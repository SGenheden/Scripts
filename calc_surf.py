# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to calculate molecular surfaces and roughness thereof

I needs the MSMS software installed
"""

import subprocess
import tempfile
import os

import numpy as np

# THIS NEED TO BE UPDATED
MSMS_PATH = "/home/sg6e12/Programs/msms.x86_64Linux2.2.6.1"

bornradii = {"H":1.2,"C":1.7,"N":1.55,"O":1.5,"NA":1.2,"P":1.8,"S":1.8,"CL":1.7,"FE":1.3}

def calc_surf(xyzrname,natoms,probe) :
  """
  Calls MSMS to calculate surfaces

  Parameters
  ----------
  xyzrname : string
    name of file containing Cartesian coordinates and radii
  natoms : int
    the number of atoms in the system
  probe : float
    the probe radius for the surface calculation

  Returns
  -------
  numpy array
    the SES and SAS for each atom
  """
  msmsfile,msmsname = tempfile.mkstemp()
  surfname = "calcsurf_surffile.area"

  command = "%s -if %s -prob %.3f -af %s"%(MSMS_PATH,xyzrname,probe,surfname)
  ret_code = subprocess.call ( command, shell = True , stdout=msmsfile, stderr=msmsfile)

  surf = np.zeros([natoms,2])
  with open(surfname,"r") as f :
    dummy = f.readline()
    for i in range(natoms) :
      line = f.readline()
      if not line : break
      surf[i,:] = map(float,line.strip().split()[1:])

  os.remove(msmsname)
  os.remove(surfname)
  return surf

def write_xyzr(xyz,radius) :
  """
  Write Cartesian coordinates + radii to file

  The file is a temporary file that should be deleted by
  the calling function

  Parameters
  ----------
  xyz : Numpy array
    the Cartesian coordinates
  radius : Numpy array
    the radii

  Returns
  -------
  string
    the name of the temporary file
  """
  xyzrfile,xyzrname = tempfile.mkstemp()
  with os.fdopen(xyzrfile,"w") as f :
    for coord,rad in zip(xyz,radius) :
      f.write("%.3f %.3f %.3f %.3f\n"%(coord[0],coord[1],coord[2],rad))
  return xyzrname

class Selection :
  """
  Class to store a residue selection
  that are used to calculate surfaces

  Attributes
  ----------
  residues : list of Residue objects
    the residues
  atoms : list of Atom object
    the atoms of all residues
  """
  def __init__(self,residues,xyzrname) :
    self.residues = residues
    self.atoms = []
    for residue in residues :
      self.atoms.extend(residue.atoms)
    self.xyzrname = xyzrname

    # Assign which atoms to sum up for each residue
    self._restake = np.zeros([len(residues)-8,len(self.atoms)],bool)
    for i,res in enumerate(self.residues[4:-4]) :
      ii = i + 4
      # Take atoms N-4 and N-3
      resib1 = max(0,ii-4)
      resib2 = max(0,ii-3)
      ai1 = self.residues[resib1].atoms[0].idx-self.atoms[0].idx
      ai2 = self.residues[resib2].atoms[-1].idx-self.atoms[0].idx
      self._restake[i,ai1:ai2+1] = True

      # Take atoms N+3 and N+4
      resif1 = min(len(self.residues)-1,ii+3)
      resif2 = min(len(self.residues)-1,ii+4)
      ai1 = self.residues[resif1].atoms[0].idx-self.atoms[0].idx
      ai2 = self.residues[resif2].atoms[-1].idx-self.atoms[0].idx
      self._restake[i,ai1:ai2+1] = True

      # Take N
      ai1 = self.residues[i].atoms[0].idx-self.atoms[0].idx
      ai2 = self.residues[i].atoms[-1].idx-self.atoms[0].idx
      self._restake[i,ai1:ai2+1] = True
    # Here we will store the residue contribution to the area
    self._patchareas = []
    self._resareas = []
    # Probe radii
    self._probes = []
  def calc_surf(self,probe) :
    areas = calc_surf(self.xyzrname,len(self.atoms),probe)
    patcharea = np.zeros([len(self._restake),2])
    for i,take in enumerate(self._restake) :
      patcharea[i,:] = areas[take,:].sum(axis=0)
    self._patchareas.append(patcharea)
    resarea = np.zeros([len(self.residues),2])
    for i,res in enumerate(self.residues) :
      ai1 = res.atoms[0].idx-self.atoms[0].idx
      ai2 = res.atoms[-1].idx-self.atoms[0].idx
      resarea[i,:] = areas[ai1:ai2+1,:].sum(axis=0)
    self._resareas.append(resarea)
    self._probes.append(probe)
  def patch_areas(self) :
    return np.array(self._patchareas)
  def res_areas(self) :
    return np.array(self._resareas)
  def fractal(self) :
    patch_areas = self.patch_areas()
    d = np.zeros([2,len(self._restake)])
    lnprobes = np.log(np.array(self._probes))
    for i,(res,take) in enumerate(zip(self.residues[4:-4],self._restake)) :
      d[0,i] = res.idx
      d[1,i] = 2.0 - np.polyfit(lnprobes,np.log(patch_areas[:,i,0]),1)[0]
    return d

if __name__ == '__main__' :

  import sys
  import argparse
  import matplotlib.pylab as plt

  from Pdb import pdb

  # Command-line input
  parser = argparse.ArgumentParser(description="Calculating surfaces")
  parser.add_argument('-f','--file',help="pdb file to analyze",default="")
  parser.add_argument('-p','--probes',nargs="+",type=float,help="probe size(s)",default=[1.4])
  parser.add_argument('-s','--selections',nargs="+",help="sub selection to make",default=[])
  args = parser.parse_args()

  pdbfile = pdb.PDBFile(filename=args.file)

  # Assign Born radii to each atom
  radius = np.zeros(pdbfile.xyz.shape[0])
  for i,atom in enumerate(pdbfile.atoms) :
    radius[i] = bornradii[atom.element().upper()]

  # Write total system and selections to individual PDB-files
  xyzrnames = [None]*(len(args.selections)+1)
  xyzrnames[-1] = write_xyzr(pdbfile.xyz,radius)
  selections = [None]*len(args.selections)
  for i,sel in enumerate(args.selections) :
    seli = map(int,sel.split("-"))
    aidx1 = pdbfile.residues[seli[0]-1].atoms[0].idx
    aidx2 = pdbfile.residues[seli[1]-1].atoms[-1].idx
    xyzrnames[i] = write_xyzr(pdbfile.xyz[aidx1:aidx2+1,:],radius[aidx1:aidx2+1])
    selections[i] = Selection(pdbfile.residues[seli[0]:seli[1]+1],xyzrnames[i])

  calcfractal = len(args.probes)>1
  if not calcfractal : print "SESA SASA"
  for i,probe in enumerate(args.probes) :
    areas = calc_surf(xyzrnames[-1],len(pdbfile.atoms),probe)
    total = areas.sum(axis=0)
    print "%.3f %.3f"%(total[0],total[1])
    for sel in selections : sel.calc_surf(probe)

  print ""
  if calcfractal :
    print "Fractal:"
    means = []
    for j,sel in enumerate(selections) :
      fractal = sel.fractal()
      plt.subplot(2,int(np.ceil(len(args.selections)/2.0)),j+1)
      plt.plot(fractal[0,:],fractal[1,:])
      print "%d-%d : %.4f +- %.4f"%(sel.residues[0].idx,sel.residues[-1].idx,fractal[1,:].mean(),fractal[1,:].std())
      means.append(fractal[1,:].mean())
    means = np.array(means)
    print "On average: %.4f +- %.4f"%(means.mean(),means.std())
    plt.tight_layout()
    #plt.show()
  else :
    print "Exposure:"
    means = []
    for j,sel in enumerate(selections) :
      free_area = sel.res_areas()[0,:,0]
      rsel = free_area>0.0
      real_area = np.zeros(len(sel.residues))
      for i,res in enumerate(sel.residues) :
        ai1 = res.atoms[0].idx-sel.atoms[0].idx
        ai2 = res.atoms[-1].idx-sel.atoms[0].idx
        real_area[i] = areas[ai1:ai2+1,0].sum()
      fraction = real_area[rsel] / free_area[rsel]
      print "%d-%d : %.4f +- %.4f"%(sel.residues[0].idx,sel.residues[-1].idx,fraction.mean(),fraction.std())
      means.append(fraction.mean())
      plt.subplot(2,int(np.ceil(len(args.selections)/2.0)),j+1)
      x = np.arange(sel.residues[0].idx,sel.residues[-1].idx+1)
      x = x[rsel]
      plt.plot(x,fraction)#
    means = np.array(means)
    print "On average: %.4f +- %.4f"%(means.mean(),means.std())
    plt.tight_layout()
    #plt.show()

  # Remove temporary files
  for xyzrname in xyzrnames :
    os.remove(xyzrname)
