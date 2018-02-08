# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to draw roughness (fractal) pies on structures

No arguments are necessary, all structures are taken from standard locations
"""

import argparse
import os
import sys

import numpy as np

import matplotlib
if not "DISPLAY" in os.environ or os.environ["DISPLAY"] == "" :
  matplotlib.use('Agg')
import matplotlib.pylab as plt
import matplotlib.colors as colors
import matplotlib.patches as patches

import  gpcr_lib
# Import the calc_surf program
thispath = os.path.dirname(os.path.abspath(__file__))
oneup = os.path.split(thispath)[0]
sys.path.insert(0,oneup)
import calc_surf

class StructurePartition :
  """
  Class to store a structure partition on pies
  to perform fractal calculations

  Attributes
  ----------
  aastruct : PDBFile
    the atomistic structure used in the surface calculations
  edges : numpy array
    the radial edges of the partition
  fractal_low : numpy array
    the average fractal for each partition in the lower leaflet
    is None unless calc_fractal() has been called
  fractal_upp : numpy array
    the average fractal for each partition in the upper leaflet
    is None unless calc_fractal() has been called
  lowsel : numpy array
    indicates if residue is on lower leaflet
  minval : float
    the lowest fractal
  maxval : float
    the highest fractal
  npies : float
   the number of radial partitions, i.e. pies
  probes : numpy array
    the probe radius
    is None if calc_surface() has not been called
  selection : list of Selection objects
    selection for each helix
  xray : XrayDensity object
    the structure of the X-ray
  xyzrnames : list of string
    filename for temporary structure names
  """
  def __init__(self, xray, aastruct, npies) :
    self.xray = xray
    self.npies = npies
    self.aastruct = aastruct

    # Determine in which pie the residues are in
    allcent = xray.pdbfile.xyz.mean(axis=0)
    centres = np.array([res.collect("centerofmass")-allcent for res in xray.pdbfile.residues])
    self.edges = np.linspace(-180,180,self.npies+1,endpoint=True)
    ang = np.arctan2(centres[:,0],centres[:,1])*180.0/np.pi
    self.partition = np.digitize(ang,self.edges)

    # Setup a calc_surf.Selection for each helix
    self.xyzrnames = [None]*len(xray.template.rhelices)
    self.selections = [None]*len(xray.template.rhelices)
    radii = np.asarray([calc_surf.bornradii[atom.element().upper()] for atom in aastruct.atoms])
    for i,h in enumerate(xray.template.rhelices) :
      aidx1 = aastruct.residues[h[0]-1].atoms[0].idx
      aidx2 = aastruct.residues[h[1]-1].atoms[-1].idx
      self.xyzrnames[i] = calc_surf.write_xyzr(aastruct.xyz[aidx1:aidx2+1,:],radii[aidx1:aidx2+1])
      self.selections[i] = calc_surf.Selection(aastruct.residues[h[0]:h[1]+1],self.xyzrnames[i])

    # This select residues in lower leaflet
    self.lowsel = centres[:,2]+allcent[2] < xray.box[2] / 2.0

    # Initialise arrays to None
    self.probes = None
    self.fractal_low = None
    self.fractal_upp = None

  def clean_up(self) :
    """
    Removes the temporary structure files from disc
    """
    for xyzrname in self.xyzrnames :
      os.remove(xyzrname)

  def calc_fractal(self) :
    """
    Calculates the fractal for each partion, by averaging over helix selections
    """
    if self.probes is None : return

    self.fractal_low = np.zeros(self.npies)
    self.fractal_upp = np.zeros(self.npies)
    ncount_low = np.zeros(self.npies)
    ncount_upp = np.zeros(self.npies)
    for sel in self.selections :
      fsel =  sel.fractal()
      for f in fsel.T :
        res = int(f[0])
        part = self.partition[res] - 1
        if self.lowsel[res] :
          self.fractal_low[part] += f[1]
          ncount_low[part] += 1.0
        else :
          self.fractal_upp[part] += f[1]
          ncount_upp[part] += 1.0
    self.fractal_low = self.fractal_low / ncount_low
    self.fractal_upp = self.fractal_upp / ncount_upp
    self.minval = min(self.fractal_low.min(),self.fractal_upp.min())
    self.maxval = max(self.fractal_low.max(),self.fractal_upp.max())

  def calc_surf(self,probes) :
    """
    Calculates the surface of the helix selections for different probe radii
    """
    self.probes = probes
    for probe in probes :
      for sel in self.selections :
        sel.calc_surf(probe)

  def plot(self, axes, labels, restocolor) :
    """
    Plot the fractal on a pie chart

    Parameters
    axes : tuple of Axis object
      the axis to draw the lower and upper leaflet pie chart
    labels : tuple of strings
      the labels to draw next to each axis
    restocolor : list of int
      residues to color
    """

    def draw_pies(axis,fractal,cent,rad,reverseY) :

      c = plt.Circle(cent.T,rad,ec='k',fc=None,fill=False)
      axis.add_patch(c)
      fractal2 = (fractal - self.minval) / (self.maxval - self.minval)
      for val,e1,e2 in zip(fractal2,self.edges[:-1],self.edges[1:]) :
        w = patches.Wedge(cent.T,rad,e1,e2,ec='k',fc='k',alpha=0.1)
        axis.add_patch(w)
        if reverseY :
          x = -rad*np.cos(e1*np.pi/180.0)
          y = rad*np.sin(e1*np.pi/180.0)
          ee2 = 180.0*np.arctan2(y,x)/np.pi
          x = -rad*np.cos(e2*np.pi/180.0)
          y = rad*np.sin(e2*np.pi/180.0)
          ee1 =  180.0*np.arctan2(y,x)/np.pi
          w = patches.Wedge(cent.T,rad,ee1,ee2,ec=plt.cm.RdYlBu_r(val),fc=plt.cm.RdYlBu_r(val),width=5)
        else :
          w = patches.Wedge(cent.T,rad,e1,e2,ec=plt.cm.RdYlBu_r(val),fc=plt.cm.RdYlBu_r(val),width=5)
        axis.add_patch(w)

    if self.fractal_low is None : return

    gpcr_lib.plot_density_xray(axes[0],0,"",0,0,self.xray,"low","Intra.",number=None,plotn=False,drawchol=False, specialres=restocolor)
    gpcr_lib.plot_density_xray(axes[1],0,"",0,0,self.xray,"upp","Extra.",number=None,plotn=False,drawchol=False, specialres=restocolor)

    rad = 30.0
    cent = np.array([0.0,0.0])
    draw_pies(axes[0],self.fractal_low,cent,rad,False)
    draw_pies(axes[1],self.fractal_upp,cent,rad,True)

    for a,l in zip(axes,labels) :
      a.text(-40,38,l)
      a.set_xticklabels([])
      a.set_yticklabels([])

  def print_helixroughness(self) :
    all = []
    for i, sel in enumerate(self.selections,1):
      fsel = sel.fractal()
      av = fsel[1,:].mean()
      print "\tH%d\t%.5f\t%.5f"%(i, av, fsel[1,:].std())
      all.append(av)
    all = np.asarray(all)
    print "\tOverall\t%.5f\t%.5f"%(all.mean(),all.std()/np.sqrt(all.shape[0]))

if __name__ == '__main__' :

  # Command-line input
  parser = argparse.ArgumentParser(description="Plotting fractal pies")
  parser.add_argument('-f','--folder', help="the  folder with the residue contacts")
  parser.add_argument('-n','--npies',type=int,help="the number of pies",default=12)
  parser.add_argument('-p','--probes',nargs="+",type=float,help="the probe sizes",default=[1.4,1.8,2.2,2.6,3.0])
  args = parser.parse_args()

  mols = "b2 b2_a a2a a2a_a".split()
  numbers = "A) B) C) D) E) F) G) H)".split()
  fig = plt.figure(1,figsize=(8,12))

  # Setup and calculate the partition for each molecule
  parts = [None]*len(mols)
  for i,mol in enumerate(mols) :
    xray, aastruct = gpcr_lib.load_xray(mol, loadsigma=True, loadaa=True)
    parts[i] = StructurePartition(xray, aastruct, args.npies)
    parts[i].calc_surf(args.probes)
    parts[i].calc_fractal()
    parts[i].clean_up()

  # Find the lowest and maxium fractal among the different molecules
  minval = min(2.0,[p.minval for p in parts])
  maxval = np.around(max([p.maxval for p in parts]),1)
  # Plot each of the pie charts
  for i, (part, mol) in enumerate(zip(parts,mols)) :
    part.minval = minval
    part.maxval = maxval
    a1 = fig.add_subplot(len(mols),2,i*2+1)
    a2 = fig.add_subplot(len(mols),2,i*2+2)
    restocolor = gpcr_lib.read_rescontacts(args.folder, mol)
    part.plot([a1,a2], numbers[(i*2):(i+1)*2], restocolor)


    # This adds a colormap
    if i == 0 :
      fig_dummy = plt.figure(2)
      im = np.outer(np.arange(part.minval,part.maxval,0.01),np.ones(10))
      a = fig_dummy.add_subplot(1,1,1,aspect="equal")
      im = a.imshow(im,aspect=0.1,cmap=plt.cm.RdYlBu_r,origin='lower',extent=(0,1,part.minval,part.maxval))
      a.get_xaxis().set_visible(False)
      gpcr_lib.draw_colormap(fig,im,text="   Fractal", unittxt="")

  fig.savefig("roughness_anal.png",format="png")

  # Print average helix roughness
  for part, mol in zip(parts, mols):
      print "Average helix roughness for %s"%mol
      part.print_helixroughness()
