# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Classes and routines to read, write and manipulate LAMMPS files

Note that only a sub-set of the versatile datafile and include
file can be read, written and manipulated.

The selection has been made on a need-basis to handle ELBA force field
and dual resolution techniques.
"""

import sys
import os
import copy
from operator import attrgetter
from ConfigParser import SafeConfigParser

import numpy as np

from sgenlib import geo

def get_filename(filename) :
  """
  Returns a filename in the directory of this python file
  """
  thispath = os.path.dirname(os.path.abspath(__file__))
  return os.path.join(thispath,filename)

def parse_molecules(data,order=False) :
  """
  Parse a datafile into a dictionary of molecules

  Parameters
  ----------
  data : Datafile object
    the datafile
  order : boolean, optional
    if the order of the molecules should be returned

  Returns
  -------
  dictionary with list of Atom objects
    the molecules
  """
  molecules = {}
  if order : mollist = []
  for atom in data.atoms :
    if atom.molecule not in molecules :
      molecules[atom.molecule] = []
      if order : mollist.append(atom.molecule)
    molecules[atom.molecule].append(atom)
  if order :
    return molecules,mollist
  else :
    return molecules

def comp_pair(p1,p2) :
  """
  Compare pair parameters, and determine which
  has precedence in a sorted list
  """
  def comp_atoms(p1,p2) :
    if p1.iatom < p2.iatom :
      return -1
    elif p1.iatom > p2.iatom :
      return 1
    else :
      if p1.jatom < p2.jatom :
        return -1
      elif p1.jatom > p2.jatom :
        return 1
      else :
        return 0
  elba_func = "lj/sf/dipole/sf"
  if p1.func == elba_func  and p2.func != elba_func :
    return -1
  elif p1.func != elba_func  and p2.func == elba_func :
    return 1
  else :
    if p1.hybrid > -1 :
      if p1.hybrid < p2.hybrid :
        return 1
      elif p1.hybrid > p2.hybrid :
        return -1
      else :
        return comp_atoms(p1,p2)
    else :
      return comp_atoms(p1,p2)

class Atom :
  """
  Class to store a datafile Atom record

  Attributes
  ----------
  idx : integer
    the serial index
  atype : integer
    the atom type
  x : float
    the Cartesian x-coordinate
  y : float
    the Cartesian y-coordinate
  z : float
    the Cartesian z-coordinate
  q : float
    the charge
  molecule : integer
    the molecule index
  diameter : float
    the diameter
  density : float
    the atom density
  mux : float
    the x-component of the dipole moment
  muy : float
    the y-component of the dipole moment
  muz : float
    the z-component of the dipole moment
  xyz : Numpy array
    the Cartesian coordinates
  mu : Numpy array
    the dipole moment
  kind : string
    atom kind specifier, can be aa, cg or cg/aa
  comment : string
    a comment to append at the end of the record
  ix : integer
    the x-component of the periodic box specification
  iy : integer
    the y-component of the periodic box specification
  iy : integer
    the z-component of the periodic box specification
  fx : float
    the x-component of the force
  fy : float
    the y-component of the force
  fz : float
    the z-component of the force
  """
  def __init__(self,record=None) :
    self.idx = 0
    self.atype = 0
    self.x = 0.0
    self.y = 0.0
    self.z = 0.0
    self.q = 0.0
    self.molecule = 0
    self.diameter = None
    self.density = None
    self.mux = 0.0
    self.muy = 0.0
    self.muz = 0.0
    self.xyz = None
    self.mu = None
    self.kind = "aa"
    self.comment = ""
    self.ix = None
    self.iy = None
    self.iz = None
    if record != None : self.read(record)
  def read(self,record) :
    """
    Read a line from a Datafile
    """
    if record.find("#") > -1  :
      record,self.comment = record.split("#")
      self.comment = "# "+self.comment.strip()
    cols = record.strip().split()
    nints = 0
    try :
      test = int(cols[-1])
      nints = nints + 1
      test = int(cols[-2])
      nints = nints + 1
      test = int(cols[-3])
      nints = nints + 1
    except :
      pass
    if nints == 3 : # If ends if ix,iy,iz
       self.iz = int(cols[-1])
       self.iy = int(cols[-2])
       self.ix = int(cols[-3])
       cols = cols[:-3]

    if len(cols) == 7 :
      self.kind = "aa"
      self.set_xyz([float(cols[4]),float(cols[5]),float(cols[6])])
      self.q = float(cols[3])
      self.atype = int(cols[2])
      self.molecule = int(cols[1])
      self.diameter = None
      self.density = None
      self.set_mu([None,None,None])
    elif len(cols) == 12 :
      try :
        test = int(cols[5])
        # If the sixth item is an integer, we are dealing with a CG atom
        self.kind = "cg"
        self.set_xyz([float(cols[2]),float(cols[3]),float(cols[4])])
        self.q = float(cols[6])
        self.atype = int(cols[1])
        self.molecule = int(cols[5])
        self.diameter = float(cols[10])
        self.density = float(cols[11])
        self.set_mu([float(cols[7]),float(cols[8]),float(cols[9])])
      except :
        # If the sixth item is an integer, we are dealing with a CG/AA atom
        self.kind = "cg/aa"
        self.set_xyz([float(cols[2]),float(cols[3]),float(cols[4])])
        self.q = float(cols[7])
        self.atype = int(cols[1])
        self.molecule = int(cols[11])
        self.diameter = float(cols[5])
        self.density = float(cols[6])
        self.set_mu([float(cols[8]),float(cols[9]),float(cols[10])])
    else :
      raise Exception("Unknown atom line!")
    self.idx = int(cols[0])
  def __str__(self) :
    istr = ""
    if self.ix != None :
      istr = " %4d %4d %4d"%(self.ix,self.iy,self.iz)
    if self.kind == "aa" :
      return "%6d %4d %2d %8.5f %12.5f %12.5f %12.5f%s %s"%(self.idx,self.molecule,self.atype,self.q,self.x,self.y,self.z,istr,self.comment)
    elif self.kind == "cg" :
      return "%6d %2d %12.5f %12.5f %12.5f %4d %8.5f %12.5f %12.5f %12.5f %12.5f %12.5f%s %s"%(self.idx,self.atype,self.x,self.y,self.z,self.molecule,self.q,self.mux,self.muy,self.muz,self.diameter,self.density,istr,self.comment)
    elif self.kind == "cg/aa" :
      return "%6d %2d %12.5f %12.5f %12.5f %12.5f %12.5f %8.5f %12.5f %12.5f %12.5f %4d%s %s"%(self.idx,self.atype,self.x,self.y,self.z,self.diameter,self.density,self.q,self.mux,self.muy,self.muz,self.molecule,istr,self.comment)
  def set_xyz(self,xyz) :
    """
    Sets the Cartesian coordinates
    """
    self.x = xyz[0]
    self.y = xyz[1]
    self.z = xyz[2]
    self.xyz = np.array(xyz)
  def set_mu(self,mu) :
    """
    Sets the dipole moment
    """
    if mu[0] != None :
      self.mux = mu[0]
      self.muy = mu[1]
      self.muz = mu[2]
      self.mu = np.array(mu)
    else :
      self.mux = self.muy = self.muz = self.mu = None

class Connectivity :
  """
  Class to store a bond, angle or torsion

  Attributes
  ----------
  idx : int
    the serial index
  param : int
    the parameter index
  atoms : list of int
    the atom indices involved
  comment : string
    a comment to append at the end
  """
  def __init__(self,record=None) :
    self.idx = 0
    self.param = 0
    self.atoms = []
    self.comment = ""
    if record != None  : self.read(record)
  def offset_atoms(self,offset) :
    """
    Increase all atom indices involved in this connectivity
    """
    self.atoms = [i+offset for i in self.atoms]
  def read(self,record) :
    """
    Read a line from a Datafile
    """
    # Take care of the comment first
    if record.find("#") > -1  :
      record,self.comment = record.split("#")
      self.comment = "# "+self.comment.strip()

    cols = record.strip().split()
    self.idx = int(cols[0])
    self.param = int(cols[1])
    self.atoms = [int(i) for i in cols[2:]]
    if len(self.atoms) == 3 and self.atoms[1] == self.atoms[2] :
        self.atoms[2] = 1
        self.comment += " # modified to be compatible with latest LAMMPS version"

  def __str__(self) :
    return "%6d %2d %s %s"%(self.idx,self.param," ".join(["%6d"%i for i in self.atoms]),self.comment)
  def __cmp__(self,other) :
    for a1,a2 in zip(self.atoms,other.atoms) :
      if a1 < a2 :
        return -1
      elif a1 > a2 :
        return 1
    return 0
  def __getitem__(self,idx) :
    if idx == 0 :
      return self.param
    else :
      return self.atoms[idx-1]

class AtomType :
  """
  Class to store an atom type

  Attributes
  ----------
  idx : int
    the serial index
  mass : float
    the mass
  func : string
    the pair function
  comment : string
    a comment to append to the end
  epsilon : float
    an LJ parameter
  sigma : float
    an LJ parameter
  """
  def __init__(self) :
    self.idx = -1
    self.mass = 0.0
    self.func = ""
    self.comment = ""
    self.epsilon = None
    self.sigma = None
  def read_mass(self,record) :
    """
    Read a mass record from either a Datafile or an Includefile
    """
    cols = record.strip().split()
    off = 0
    if cols[0] == "mass" :
      off = 1
    idx = int(cols[0+off])
    if self.idx == -1 or idx == self.idx :
      self.idx = idx
      self.mass = mass = float(cols[1+off])
      if len(cols) > (2+off) :
        self.comment = " ".join(cols[(2+off):])
    else :
      print "Warning: Trying to set mass of an illegal atom type!"
  def read_lj(self,record) :
    """
    Read LJ coefficients from a Datafile
    """
    cols = record.strip().split()
    idx = int(cols[0])
    if self.idx == -1 or idx == self.idx :
      self.idx = idx
      self.epsilon = float(cols[1])
      self.sigma = float(cols[2])
      if len(cols) > 2 :
        self.comment = " ".join(cols[3:])
    else :
      print "Warning: Trying to set mass of an illegal atom type!"

class PairParam :
  """
  Class to store pair parameters for an Includefile

  Attributes
  ----------
  iatom : int
    atom type index
  jatom : int
    atom type index
  hybrid : int
    hybrid flag
  func : string
    pair function
  sigma : float
    LJ parameter
  epsilon : float
    LJ parameter
  sigma14 : float
    LJ parameter
  epsilon14 : float
    LJ parameter
  scale : float
    scaling of Coulomb
  comment : string
    a comment to append
  """
  def __init__(self,record=None) :
    self.iatom = 0
    self.jatom = 0
    self.hybrid = -1
    self.func = ""
    self.sigma = 0.0
    self.epsilon = 0.0
    self.sigma14 = 0.0
    self.epsilon14 = 0.0
    self.scale = None
    self.comment = "" # Comment
    self.ljres = 5
    if record != None  : self.read(record)
  def read(self,record) :
    """
    Read pair parameters from an Includefile or Datafile
    """
    cols = record.strip().split()
    if cols[0] != "pair_coeff" : cols.insert(0,"pair_coeff") # Needed for datafile

    self.iatom = int(cols[1])
    self.jatom = int(cols[2])
    try :
      test = float(cols[3])
      off = 0
    except :
      self.func = cols[3]
      off = 1
    try :
      test = int(cols[3+off])
      self.hybrid = test
      off = off + 1
    except :
      pass
    self.epsilon = float(cols[3+off])
    self.sigma = float(cols[4+off])
    try :
      test = float(cols[5+off])
      if self.func == "lj/sf/dipole/sf" :
        self.scale = float(cols[6+off])
        off = off + 2
      else :
        self.epsilon14 = float(cols[5+off])
        self.sigma14 = float(cols[6+off])
        off = off + 2
    except :
      pass
    if 4+off+1 < len(cols) :
      self.comment = " ".join(cols[4+off+1:])
  def __str__(self) :
    str14 = ""
    if self.func == "lj/sf/dipole/sf" and self.scale is not None:
      str14 = " scale %15.3f"%self.scale
    else :
      if self.sigma14 != 0.0 :
        str14 = " %15.10f %15.10f"%(self.epsilon14,self.sigma14)

    ljfrmstr = "%"+"15.%df"%self.ljres
    ljfrmstr = ljfrmstr + " " + ljfrmstr
    ljstr = ljfrmstr%(self.epsilon,self.sigma)

    if self.hybrid > -1 :
      return "%2d %2d %-20s %2d %s%s %s"%(self.iatom,self.jatom,self.func,self.hybrid,ljstr,str14,self.comment)
    else :
      return "%2d %2d %-20s %2s %s%s %s"%(self.iatom,self.jatom,self.func,"",ljstr,str14,self.comment)

class ConnectivityParam :
  """
  Class to store a connectivity parameter

  Attributes
  ----------
  idx : int
    the serial index
  func : string
    the functional form
  comment : string
    a comment to append
  params : list of string
    the parameters
  """
  def __init__(self,record=None) :
    self.idx = -1
    self.func = ""
    self.comment = ""
    self.params = []
    if record != None  : self.read(record)
  def read(self,record) :
    """
    Read a line from either a Datafile or an Includefile
    """
    cols = record.strip().split()
    try :
      # If the first line is an integer, we are reading from a data file
      test = int(cols[0])
      off = 0
    except :
      off = 1
    self.idx = int(cols[0+off])
    try :
      # If the next col can be converted to a string, the param function is not specified
      test = float(cols[1+off])
    except :
      self.func = cols[1+off]
      off = off + 1
    nparam = 0
    for c in cols[1+off:] :
      if c[0] == "#" : break
      nparam = nparam + 1
    self.params = [self.__parse_param(c) for c in cols[1+off:1+nparam+off]]
    if nparam+off+1 < len(cols) :
      self.comment = " ".join(cols[1+nparam+off:])
  def __parse_param(self,param) :
    """
    Parse a parameter to an integer or float and convert it to a formated string
    """
    try :
      test = int(param)
      return "%8d"%test
    except :
      test = float(param)
      return "%8.4f"%test
  def __str__(self) :
    return "%2d %s %s %s"%(self.idx,self.func," ".join(self.params),self.comment)

class Datafile :
  """
  Class to store a LAMMPS Datafile

  Attributes
  ----------
  atoms : list of Atom objects
    the atoms
  bonds : list of Connectivity objects
    the bonds
  angles : list of Connectivity objects
    the angles
  dihedrals : list of Connectivity objects
    the dihedrals
  impropers : list of Connectivity objects
    the impropers
  atomtypes : list of AtomType objects
    the atom types
  pairtypes : list of PairParam objects
    the pair params, read from PairIJ coeff
  bondtypes : list of ConnectivityParam
    the bond types
  angletypes : list of ConnectivityParam
    the angle types
  dihedraltypes : list of ConnectivityParam
    the dihedral types
  impropertypes : list of ConnectivityParam
    the improper types
  box : Numpy array
    the box coordinates
  title : string
    the title of the file
  """
  # Constructor
  def __init__(self,filename=None) :
    self.atoms = []
    self.bonds = []
    self.angles = []
    self.dihedrals = []
    self.impropers = []
    self.atomtypes = []
    self.pairtypes = []
    self.bondtypes = []
    self.angletypes = []
    self.dihedraltypes = []
    self.impropertypes = []
    self.box = np.zeros(6)
    self.title = ""
    if filename != None : self.read(filename)
  def extend(self,copy) :
    """
    Extend atoms and connectivity from another Datafile
    """
    self.title = self.title + "//" + copy.title
    # Extend atoms
    nmol = 0
    natoms = len(self.atoms)
    natomtypes = len(self.atomtypes)
    self.atomtypes.extend([None]*len(copy.atomtypes))
    for a in self.atoms : nmol = max(nmol,a.molecule)
    has_aa = self.atoms[0].kind == "aa"
    for a in copy.atoms :
      self.atoms.append(a)
      self.atoms[-1].idx = self.atoms[-1].idx + natoms
      self.atoms[-1].molecule = self.atoms[-1].molecule + nmol + 1
      self.atoms[-1].atype = self.atoms[-1].atype + natomtypes
      if not has_aa and self.atoms[-1].kind == "aa" :
        self.atoms[-1].kind = self.atoms[0].kind
        self.atoms[-1].diameter = 0.0
        self.atoms[-1].density = copy.atomtypes[self.atoms[-1].atype-1-natomtypes].mass
        self.atoms[-1].set_mu([0.0,0.0,0.0])
    # Extend connectivity
    for b in copy.bonds :
      self.bonds.append(b)
      self.bonds[-1].idx = len(self.bonds)
      self.bonds[-1].param = self.bonds[-1].param + len(self.bondtypes)
      self.bonds[-1].offset_atoms(natoms)
    self.bondtypes.extend([None]*len(copy.bondtypes))
    for a in copy.angles :
      self.angles.append(a)
      self.angles[-1].idx = len(self.angles)
      self.angles[-1].param = self.angles[-1].param + len(self.angletypes)
      self.angles[-1].offset_atoms(natoms)
    self.angletypes.extend([None]*len(copy.angletypes))
    for d in copy.dihedrals :
      self.dihedrals.append(d)
      self.dihedrals[-1].idx = len(self.dihedrals)
      self.dihedrals[-1].param = self.dihedrals[-1].param + len(self.dihedraltypes)
      self.dihedrals[-1].offset_atoms(natoms)
    self.dihedraltypes.extend([None]*len(copy.dihedraltypes))

  def read(self,filename) :
    """
    Read a file from disc
    """
    with open(filename,"r") as f :
      self.title = f.readline().strip()
      line = f.readline()
      while line :
        # Initialize empty arrays for atoms, connectivity and parameters
        if line.find(" atoms") > -1 : self.atoms = [None]*int(line.strip().split()[0])
        elif line.find(" bonds") > -1 : self.bonds = [None]*int(line.strip().split()[0])
        elif line.find(" angles") > -1 : self.angles = [None]*int(line.strip().split()[0])
        elif line.find(" dihedrals") > -1 : self.dihedrals = [None]*int(line.strip().split()[0])
        elif line.find(" impropers") > -1 : self.impropers = [None]*int(line.strip().split()[0])

        elif line.find(" atom types") > -1 : self.atomtypes = [None]*int(line.strip().split()[0])
        elif line.find(" bond types") > -1 : self.bondtypes = [None]*int(line.strip().split()[0])
        elif line.find(" angle types") > -1 : self.angletypes = [None]*int(line.strip().split()[0])
        elif line.find(" dihedral types") > -1 : self.dihedraltypes = [None]*int(line.strip().split()[0])

        # Read the box
        elif line.find("xlo") > -1 :
          cols = line.strip().split()
          self.box[0] = float(cols[0])
          self.box[3] = float(cols[1])
        elif line.find("ylo") > -1 :
          cols = line.strip().split()
          self.box[1] = float(cols[0])
          self.box[4] = float(cols[1])
        elif line.find("zlo") > -1 :
          cols = line.strip().split()
          self.box[2] = float(cols[0])
          self.box[5] = float(cols[1])

        # Read force field parameters
        elif line.find("Masses") > -1 :
          line = f.readline() # Dummy line
          for i in range(len(self.atomtypes)) :
            if self.atomtypes[i] == None :
              self.atomtypes[i] = AtomType()
            self.atomtypes[i].read_mass(f.readline())

        elif line.find("Pair Coeffs") > -1 :
          line = f.readline() # Dummy line
          for i in range(len(self.atomtypes)) :
            if self.atomtypes[i] == None :
              self.atomtypes[i] = AtomType()
            self.atomtypes[i].read_lj(f.readline())

        elif line.find("PairIJ Coeffs") > -1 :
          line = f.readline() # Dummy line
          line = f.readline()
          while len(line) > 4 and line.strip()[0].isdigit() :
            self.pairtypes.append(PairParam(record=line))
            line = f.readline()
          continue

        elif line.find("Bond Coeffs") > -1 :
          line = f.readline() # Dummy line
          for i in range(len(self.bondtypes)) :
            self.bondtypes[i] = ConnectivityParam(record=f.readline())

        elif line.find("Angle Coeffs") > -1 :
          line = f.readline() # Dummy line
          for i in range(len(self.angletypes)) :
            self.angletypes[i] = ConnectivityParam(record=f.readline())

        elif line.find("Dihedral Coeffs") > -1 :
          line = f.readline() # Dummy line
          for i in range(len(self.dihedraltypes)) :
            self.dihedraltypes[i] = ConnectivityParam(record=f.readline())

        # Read atoms and connectivity
        elif line.find("Atoms") > -1 :
          line = f.readline() # Dummy line
          for i in range(len(self.atoms)) :
            self.atoms[i] = Atom(record=f.readline())
          self.atoms.sort(key=attrgetter("idx"))

        elif line.find("Bonds") > -1 :
          line = f.readline() # Dummy line
          for i in range(len(self.bonds)) :
            self.bonds[i] = Connectivity(record=f.readline())

        elif line.find("Angles") > -1 :
          line = f.readline() # Dummy line
          for i in range(len(self.angles)) :
            self.angles[i] = Connectivity(record=f.readline())

        elif line.find("Dihedrals") > -1 :
          line = f.readline() # Dummy line
          for i in range(len(self.dihedrals)) :
            self.dihedrals[i] = Connectivity(record=f.readline())

        line = f.readline()

  def single_lj(self) :
    """
    Return a list of AtomType objects for eahc unique atom type among the PairParam objects
    """
    if len(self.pairtypes) == 0 :
      return self.atomtypes

    paramlist = []
    for p in self.pairtypes :
      if p.iatom != p.jatom : continue
      paramlist.append(AtomType())
      paramlist[-1].func = p.func
      paramlist[-1].comment = p.comment
      paramlist[-1].sigma = p.sigma
      paramlist[-1].epsilon = p.epsilon
      paramlist[-1].idx = p.iatom
      for t in self.atomtypes :
        if t.idx == p.iatom :
          paramlist[-1].mass = t.mass
          break
    return paramlist

  def sort(self,reorder=True) :
    """
    Sort Atoms, Bonds, Angles and Dihedrals

    Parameters
    ----------
    reorder : boolean, optional
      if true, sorts connectivity based on atom indices
    """
    self.atoms.sort(key=attrgetter("idx"))

    for list in [self.bonds,self.angles,self.dihedrals] :
      if reorder :
        list.sort()
        for i,con in enumerate(list,1) :
          con.idx = i
      else :
        list.sort(key=attrgetter("idx"))

  def write(self,filename,writeparams=False) :
    """
    Write the atoms and connectivity to disc

    It does not write force field parameters
    """
    with open(filename,"w") as f :
      f.write("%s\n"%self.title)
      f.write("\n")
      f.write("%d atoms\n"%len(self.atoms))
      f.write("%d bonds\n"%len(self.bonds))
      f.write("%d angles\n"%len(self.angles))
      f.write("%d dihedrals\n"%len(self.dihedrals))
      f.write("%d impropers\n"%len(self.impropers))
      f.write("\n")
      f.write("%d atom types\n"%len(self.atomtypes))
      f.write("%d bond types\n"%len(self.bondtypes))
      f.write("%d angle types\n"%len(self.angletypes))
      f.write("%d dihedral types\n"%len(self.dihedraltypes))
      f.write("%d improper types\n"%len(self.impropertypes))
      f.write("\n")
      f.write("%15.5f %15.5f xlo xhi\n"%(self.box[0],self.box[3]))
      f.write("%15.5f %15.5f ylo yhi\n"%(self.box[1],self.box[4]))
      f.write("%15.5f %15.5f zlo zhi\n"%(self.box[2],self.box[5]))
      f.write("\n")

      if writeparams :
        if len(self.atomtypes) > 0 and self.atomtypes[0] is not None :
          f.write("Masses\n\n")
          for a in self.atomtypes : f.write("%2d %10.5f %s\n"%(a.idx,a.mass,a.comment))
          f.write("\n")
          if self.atomtypes[0].sigma is not None :
            f.write("Pair Coeffs\n\n")
            for a in self.atomtypes : f.write("%2d %s %10.5f %10.5f %s\n"%(a.idx,a.func,a.epsilon,a.sigma,a.comment))
            f.write("\n")
          if len(self.pairtypes) > 0 :
            f.write("PairIJ Coeffs\n\n")
            for p in self.pairtypes : f.write("%s\n"%p.__str__())
            f.write("\n")
        if len(self.bondtypes) > 0 and self.bondtypes[0] is not None :
          f.write("Bond Coeffs\n\n")
          for c in self.bondtypes : f.write("%s\n"%c.__str__())
          f.write("\n")
        if len(self.angletypes) > 0 and self.angletypes[0] is not None :
          f.write("Angle Coeffs\n\n")
          for c in self.angletypes : f.write("%s\n"%c.__str__())
          f.write("\n")
        if len(self.dihedraltypes) > 0 and self.dihedraltypes[0] is not None :
          f.write("Dihedral Coeffs\n\n")
          for c in self.dihedraltypes : f.write("%s\n"%c.__str__())
          f.write("\n")

      f.write("Atoms\n\n")
      for atm in self.atoms :
        f.write("%s\n"%atm.__str__())
      if len(self.bonds) > 0 : f.write("\nBonds\n\n")
      for bnd in self.bonds :
        f.write("%s\n"%bnd.__str__())
      if len(self.angles) > 0 : f.write("\nAngles\n\n")
      for ang in self.angles :
        f.write("%s\n"%ang.__str__())
      if len(self.dihedrals) > 0 : f.write("\nDihedrals\n\n")
      for dih in self.dihedrals :
        f.write("%s\n"%dih.__str__())
      if len(self.impropers) > 0 : f.write("\nImpropers\n\n")
      for imp in self.impropers :
        f.write("%s\n"%imp.__str__())
      f.write("\n")

class Includefile() :
  """
  Class to store a LAMMPS Includefile

  Attributes
  ----------
  masses : list of AtomType objects
    the masses
  pair_coeff : list of PairParam objects
    the pair coefficients
  bondparams : list of ConnectivityParam objects
    the bond parameters
  angleparams : list of ConnectivityParam objects
    the angle parameters
  dihedralparams : list of ConnectivityParam objects
    the dihedral parameters
  improperparams : list of ConnectivityParam objects
    the improper parameters
  """
  def __init__(self,filename=None) :
    self.masses = []
    self.pair_coeff = []
    self.bondparams = []
    self.angleparams = []
    self.dihedralparams = []
    self.improperparams = []
    if filename != None : self.read(filename)
  def duplicate_type(self,atype) :
    """
    Duplicate a type, creating new mass and pair_coeff entries
    """
    newmass = copy.deepcopy(self.masses[atype-1])
    newmass.idx = len(self.masses)+1
    new_paircoeff = []
    for paircoeff in self.pair_coeff :
      if atype in [paircoeff.iatom , paircoeff.jatom] :
        newpair = copy.deepcopy(paircoeff)
        if newpair.comment == "" : newpair.comment = "#"
        newpair.comment = newpair.comment + " copy of %d,%d"%(newpair.iatom,newpair.jatom)
        if paircoeff.iatom == atype :
          newpair.iatom = len(self.masses)+1
        if paircoeff.jatom == atype :
          newpair.jatom = len(self.masses)+1
        if newpair.jatom < newpair.iatom :
          temp = newpair.jatom
          newpair.jatom = newpair.iatom
          newpair.iatom = temp
        new_paircoeff.append(newpair)
        if paircoeff.iatom == paircoeff.jatom :
          newpair = copy.deepcopy(paircoeff)
          newpair.comment += " copy of %d,%d"%(newpair.iatom,newpair.jatom)
          newpair.jatom = len(self.masses)+1
          new_paircoeff.append(newpair)
    return (newmass,new_paircoeff)
  def extend(self,copy,lj_hybrid=-1,lj_func="",ang_func=None,dih_func=None) :
    """
    Extend with another Include file

    Mixes LJ parameters with Lorentz-Berthelot rules

    Parameters
    ----------
    copy : Includefile object
      the include file to append to this
    lj_hybrid : int
      the hybrid flag for the pair coefficient
    lj_func : string
      the function for the pair coefficient
    ang_func : string, optional
      the function for the angle parameters
    dih_func : string, optional
      the function for the dihedral parameters
    """
    natomtypes = 0
    for m in self.masses : natomtypes = max(natomtypes,m.idx)
    # Add new masses and offset the idx
    for m in copy.masses :
      self.masses.append(m)
      self.masses[-1].idx = natomtypes + self.masses[-1].idx
    # Add new pair coefficients
    my_lj = self.__lj_same()
    copy_lj = copy.__lj_same()
    for i_lj in my_lj :
      for j_lj in copy_lj :
        self.pair_coeff.append(PairParam())
        self.pair_coeff[-1].iatom = i_lj.idx
        self.pair_coeff[-1].jatom = j_lj.idx+natomtypes
        self.pair_coeff[-1].hybrid = lj_hybrid
        self.pair_coeff[-1].func = lj_func
        self.pair_coeff[-1].comment = "# Mixed using Lorentz-Berthelot rules"
        self.pair_coeff[-1].epsilon = np.sqrt(i_lj.epsilon*j_lj.epsilon)
        self.pair_coeff[-1].sigma = (i_lj.sigma+j_lj.sigma)/2.0
    for p in copy.pair_coeff :
      self.pair_coeff.append(p)
      self.pair_coeff[-1].iatom  = natomtypes + self.pair_coeff[-1].iatom
      self.pair_coeff[-1].jatom  = natomtypes + self.pair_coeff[-1].jatom
    # Add new connectivity params
    self.__extend_con(self.bondparams,copy.bondparams)
    self.__extend_con(self.angleparams,copy.angleparams,func=ang_func)
    self.__extend_con(self.dihedralparams,copy.dihedralparams,func=dih_func)
    self.__extend_con(self.improperparams,copy.improperparams)

  def extend_from_data(self,copy,lj_hybrid=None,lj_func=None,lj_hybrid_mix=None,lj_func_mix=None,ang_func=None,dih_func=None) :
    """
    Extend with another Datafile

    Mixes LJ parameters with Lorentz-Berthelot rules

    Parameters
    ----------
    copy : Includefile object
      the include file to append to this
    lj_hybrid : int, optional
      the hybrid flag for the pair coefficients in the Datafile
    lj_func : string , optional
      the function for the pair coefficients in the Datafile
    lj_hybrid_mix : dictionary, optional
      the hybrid flag for the pair coefficients when mixing them
    lj_func : dictionary, optional
      the function for the pair coefficients when mixing them
    ang_func : string, optional
      the function for the angle parameters
    dih_func : string, optional
      the function for the dihedral parameters
    """
    natomtypes = 0
    for m in self.masses : natomtypes = max(natomtypes,m.idx)
    copy_atomtypes = copy.single_lj()
    # Add new masses and offset the idx
    for m in copy_atomtypes :
      self.masses.append(m)
      self.masses[-1].idx = natomtypes + self.masses[-1].idx
    # Add new pair coefficients
    my_lj = self.__lj_same()
    for i_lj in my_lj :
      for j_lj in copy_atomtypes :
        self.pair_coeff.append(PairParam())
        self.pair_coeff[-1].iatom = i_lj.iatom
        self.pair_coeff[-1].jatom = j_lj.idx
        if lj_hybrid_mix == None :
          self.pair_coeff[-1].hybrid = -1
        else :
          self.pair_coeff[-1].hybrid = lj_hybrid_mix[i_lj.func]
        if lj_func_mix == None :
          self.pair_coeff[-1].func = ""
        else :
          self.pair_coeff[-1].func = lj_func_mix[i_lj.func]
        self.pair_coeff[-1].comment = "# Mixed using Lorentz-Berthelot rules"
        self.pair_coeff[-1].epsilon = np.sqrt(i_lj.epsilon*j_lj.epsilon)
        self.pair_coeff[-1].sigma = (i_lj.sigma+j_lj.sigma)/2.0
    if not copy.pairtypes :
      for i_lj in copy.atomtypes :
        for j_lj in copy.atomtypes :
          if j_lj.idx < i_lj.idx : continue
          self.pair_coeff.append(PairParam())
          #print i_lj.idx,j_lj.idx,natomtypes
          self.pair_coeff[-1].iatom = i_lj.idx
          self.pair_coeff[-1].jatom = j_lj.idx
          if lj_hybrid == None :
            self.pair_coeff[-1].hybrid = -1
          else :
            self.pair_coeff[-1].hybrid = lj_hybrid
          if lj_func == None :
            self.pair_coeff[-1].func = ""
          elif self.pair_coeff[-1] != "" :
            self.pair_coeff[-1].func = lj_func
          self.pair_coeff[-1].comment = "# From Datafile, Mixed using Lorentz-Berthelot rules"
          self.pair_coeff[-1].epsilon = np.sqrt(i_lj.epsilon*j_lj.epsilon)
          self.pair_coeff[-1].sigma = (i_lj.sigma+j_lj.sigma)/2.0
    else :
      for pt in copy.pairtypes :
        self.pair_coeff.append(pt)
        self.pair_coeff[-1].iatom = self.pair_coeff[-1].iatom + natomtypes
        self.pair_coeff[-1].jatom = self.pair_coeff[-1].jatom + natomtypes
        if lj_hybrid == None :
          self.pair_coeff[-1].hybrid = -1
        else :
          self.pair_coeff[-1].hybrid = lj_hybrid
        if lj_func == None :
          self.pair_coeff[-1].func = ""
        else :
          self.pair_coeff[-1].func = lj_func
    # Add new connectivity params
    self.__extend_con(self.bondparams,copy.bondtypes)
    self.__extend_con(self.angleparams,copy.angletypes,func=ang_func)
    self.__extend_con(self.dihedralparams,copy.dihedraltypes,func=dih_func)
    self.__extend_con(self.improperparams,copy.impropertypes)
  def read(self,filename) :
    """
    Read from disc
    """
    with open(filename,"r") as f :
      line = f.readline()
      while line :
        if line.find("mass") == 0 :
          self.masses.append(AtomType())
          self.masses[-1].read_mass(line)
        elif line.find("pair_coeff") == 0 : self.pair_coeff.append(PairParam(record=line))
        elif line.find("bond_coeff") == 0 : self.bondparams.append(ConnectivityParam(record=line))
        elif line.find("angle_coeff") == 0 : self.angleparams.append(ConnectivityParam(record=line))
        elif line.find("dihedral_coeff") == 0 : self.dihedralparams.append(ConnectivityParam(record=line))
        elif line.find("improper_coeff") == 0 : self.improperparams.append(ConnectivityParam(record=line))
        line = f.readline()
  def write(self,filename,ljres=5) :
    """
    Write to disc
    """
    f = open(filename,"w")
    f.write("# Masses\n")
    f.write("#    %2s %10s\n"%("type","mass"))
    for m in self.masses : f.write("mass %2d %10.5f %s\n"%(m.idx,m.mass,m.comment))
    f.write("\n# Lennard-Jones coefficients\n")
    f.write("#          %2s %2s %-20s %2s %15s %15s\n"%("i","j","(style)","(hy)","epsilon","sigma"))
    for pair in self.pair_coeff :
        pair.ljres = ljres
        f.write("pair_coeff %s\n"%pair.__str__())
    f.write("\n# Bond coefficients\n")
    f.write("#           type params\n")
    for bnd in self.bondparams : f.write("bond_coeff %s\n"%bnd.__str__())
    f.write("\n# Angle coefficients\n")
    f.write("#           type (func) params\n")
    for ang in self.angleparams : f.write("angle_coeff %s\n"%ang.__str__())
    f.write("\n# Dihedral coefficients\n")
    f.write("#             type params\n")
    for dih in self.dihedralparams : f.write("dihedral_coeff %s\n"%dih.__str__())
    f.write("\n# Improper coefficients\n")
    f.write("#             type params\n")
    for imp in self.improperparams : f.write("improper_coeff %s\n"%imp.__str__())
    f.close()
  def __lj_same(self) :
    """
    Return a list of PairParam where i==j
    """
    lj_same = []
    for pair in self.pair_coeff :
      if pair.iatom == pair.jatom : lj_same.append(pair)
    return lj_same
  def lj_same_dict(self) :
    """
    Return a dictionary of PairParam where i==j, where the key is iatom
    """
    lj_same = {}
    for pair in self.pair_coeff :
      if pair.iatom == pair.jatom : lj_same[pair.iatom] = pair
    return lj_same
  def __extend_con(self,con1,con2,func=None) :
    """
    Extend a list of connectivity params
    """
    ntypes = 0
    for con in con1 : ntypes = max(ntypes,con.idx)
    for con in con2 :
      con1.append(con)
      con1[-1].idx = ntypes + con1[-1].idx
      if func != None and con1[-1].func == "" : con1[-1].func = func

class DumpFile :
  """
  Class to read a LAMMPS DumpFile into a Datafile

  Attributes
  ----------
  datafile : DataFile object
    a template datafile
  box : NumpyArray
    the box coordinates
  idx : int
     the snapshot serial number
  """
  def __init__(self,datafile,filename=None,snapshot=None) :
    self.datafile = datafile
    self.box = None
    self.idx = None
    if filename != None : self.read(filename,snapshot)
  def read(self,filename,snapshot) :
    """
    Read a specific snapshot from disc
    """
    if snapshot == None or snapshot < 0 :
      f = open(filename,"r")
      line = f.readline()
      while line :
        if line.find("ITEM: TIMESTEP") == 0 :
          snapshot = int(f.readline().strip())
        line = f.readline()
      f.close()

    #print "Extracting snapshot %d \nReading: "%snapshot,
    f = open(filename,"r")
    line = f.readline()
    n = 0
    while line :
      if line.find("ITEM: TIMESTEP") == 0 :
        n = n + 1
        # Read timestep
        idx = int(f.readline().strip())
        if n % 1000 == 0 :
          print idx,
          sys.stdout.flush()
        if idx == snapshot :
          box = np.zeros(6)
          atoms = []
          keys = []
          # Read number of atoms
          f.readline()
          natom = int(f.readline().strip())
          if natom != len(self.datafile.atoms) :
            print "Incorrect dump file. Not the same number of atoms!"
            f.close()
            return
          # Read box
          f.readline()
          v1,v2 = f.readline().strip().split()
          box[0] = float(v1)
          box[3] = float(v2)
          v1,v2 = f.readline().strip().split()
          box[1] = float(v1)
          box[4] = float(v2)
          v1,v2 = f.readline().strip().split()
          box[2] = float(v1)
          box[5] = float(v2)
          # Read atom info
          line = f.readline()
          keys = line.strip().split()[2:]
          for i in range(natom) :
            atoms.append(f.readline().strip().split())
          break
      line = f.readline()
    f.close()
    #print ""

    self.datafile.box = box
    for dumpatom in atoms :
      idx = int(dumpatom[0])
      dataatom = self.datafile.atoms[idx-1]
      for key,val in zip(keys,dumpatom) :
        if key == "x" :
          dataatom.x = float(val)
        elif key == "y" :
          dataatom.y = float(val)
        elif key == "z" :
          dataatom.z = float(val)
        elif key == "mux" :
          dataatom.mux = float(val)
        elif key == "muy" :
          dataatom.muy = float(val)
        elif key == "muz" :
          dataatom.muz = float(val)
        elif key == "ix" :
          dataatom.ix = int(val)
        elif key == "iy" :
          dataatom.iy = int(val)
        elif key == "iz" :
          dataatom.iz = int(val)
        elif key == "fx" :
          dataatom.fx = float(val)
        elif key == "fy" :
          dataatom.fy = float(val)
        elif key == "fz" :
          dataatom.fz = float(val)


class ResidueConverter :
  """
  Class to convert all-atom molecule to CG

  Attributes
  ----------
  cg_names : list of strings
    the name of the CG beads
  aa_names : dictionary of list of strings
    the AA names that make up each CG bead
  creates : dictionary of strings
    the string used to create the Atom object for each CG bead
  bonds : list of integers
    the bond specification
  angles : list of integers
    the angle specification
  dihedrals : list of integers
    the dihedral specification
  name : string
    the name of the molecule
  """
  def __init__(self) :
    self.cg_names = []
    self.aa_names = {}
    self.creates = {}
    self.bonds = []
    self.angles = []
    self.dihedrals = []
    self.name = ""
  def generate_cg(self,residue,molidx,data,mapping=True) :
    """
    Generate a CG molecule

    Parameters
    ----------
    residue : Residue object
      the residue to make CG
    molidx : int
      the molecular index
    data : Datafile object
      the datafile to append the CG residue to

    Returns
    -------
    numpy array
      the coordinates of the CG beads
    """

    # Find and build up the coordinates for each CG bead
    if mapping :
      coords = [np.zeros(3) for name in self.cg_names]
      counts = {name : 0 for name in self.cg_names}
      for atom in residue.atoms :
        name = atom.name.strip().lower()
        for i,cg_name in enumerate(self.cg_names) :
          if name in self.aa_names[cg_name] :
            coords[i] = coords[i] + atom.xyz
            counts[cg_name] = counts[cg_name] + 1

      for cg_name in self.cg_names :
        if counts[cg_name] != len(self.aa_names[cg_name]) and cg_name is not "w" :
          print "Could not find all atoms for %s (%d %d)"%(cg_name,counts[cg_name],len(self.aa_names[cg_name]))

    else :
      coords = [np.zeros(3) for name in self.cg_names]
      counts = {name : 1 for name in self.cg_names}
      for i,ratom in enumerate(residue.atoms) :
        coords[i] = np.array(ratom.xyz,copy=True)

    # Create the Atom objects
    n = len(data.atoms)
    for i,name in enumerate(self.cg_names) :
      coords[i] = coords[i] / counts[name]
      atom = Atom(record=self.creates[name]%(len(data.atoms)+1,coords[i][0],coords[i][1],coords[i][2],molidx))
      r = np.sqrt(np.sum(atom.mu**2))
      if r > 0 :
        vec = geo.sphere_rand(r)
        rnew = np.sqrt(np.sum(vec**2))
        if abs(rnew - r) > 0.001 : print "Failed: %.3f %.3f %.3E"%(r,rnew,rnew-r)
        atom.set_mu(vec)
      data.atoms.append(atom)

    # Create the connectivity
    for bond in self.bonds :
      b = Connectivity(record="%d %d %d %d"%(len(data.bonds)+1,bond[0],bond[1]+n,bond[2]+n))
      data.bonds.append(b)

    for angle in self.angles :
      a = Connectivity(record="%d %d %d %d %d"%(len(data.angles)+1,angle[0],angle[1]+n,angle[2]+n,angle[3]+n))
      data.angles.append(a)

    for dihedral in self.dihedrals :
      a = Connectivity(record="%d %d %d %d %d %d"%(len(data.dihedrals)+1,dihedral[0],dihedral[1]+n,dihedral[2]+n,dihedral[3]+n,dihedral[4]+n))
      data.dihedrals.append(a)

    return coords
  def parse(self,section,parser) :
    """
    Parse a ConfigParser section to read in the conversion specification
    """
    self.cg_names = parser.get(section,"names").split()
    for name in self.cg_names :
      self.aa_names[name] = parser.get(section,"%s_names"%name).split()
      self.creates[name] = parser.get(section,"%s_create"%name)
    if "bonds" in parser.options(section) :
      for bond in parser.get(section,"bonds").split() :
        self.bonds.append(map(int,bond.split(",")))
    if "angles" in parser.options(section) :
      for angle in parser.get(section,"angles").split() :
        self.angles.append(map(int,angle.split(",")))
    if "dihedrals" in parser.options(section) :
      for dihedral in parser.get(section,"dihedrals").split() :
        self.dihedrals.append(map(int,dihedral.split(",")))
    self.name = section
  def __str__(self) :
    return " ".join(self.cg_names)

class Aa2Cg :
  """
  Class to store and read a dictionary of conversion specifications

  Attributes
  ----------
  residues : list of ResidueConverter objects
    the conversion specifications
  """
  def __init__(self,filename=None) :
    self.residues = []
    if filename is not None : self.read(filename)
  def read(self,filename) :
    """
    Read the file as a configuration file
    """

    parser = SafeConfigParser()
    parser.read(filename)

    # Then parse it into ResidueConvert objects
    for section in parser.sections() :
      self.residues.append(ResidueConverter())
      self.residues[-1].parse(section,parser)

#
# Legacy code, not used!!!
#

def combine_datafiles(data1,data2,box) :

  def merge_topol(tag) :

    if data1["%ss"%tag] > 0 :
      atom_offset = data1["atoms"]
      param_offset = data1["%stypes"%tag]
      data["%sdata"%tag][:data1["%ss"%tag],:] = data1["%sdata"%tag]
      if data2["%ss"%tag]  > 0 :
        data["%sdata"%tag][data1["%ss"%tag]:,0] = data2["%sdata"%tag][:,0]+param_offset
        data["%sdata"%tag][data1["%ss"%tag]:,1:] = data2["%sdata"%tag][:,1:]+atom_offset
    else :
      atom_offset = data1["atoms"]
      data["%sdata"%tag][:,0] = data2["%sdata"%tag][:,0]
      data["%sdata"%tag][:,1:] = data2["%sdata"%tag][:,1:]+atom_offset


  data = {}
  data["title"] = data1["title"] + " // " + data2["title"]
  data["atoms"] = data1["atoms"] + data2["atoms"]
  data["bonds"] = data1["bonds"] + data2["bonds"]
  data["angles"] = data1["angles"] + data2["angles"]
  data["dihedrals"] = data1["dihedrals"] + data2["dihedrals"]
  data["impropers"] = data1["impropers"] + data2["impropers"]
  data["atomtypes"] = data1["atomtypes"] + data2["atomtypes"]
  data["bondtypes"] = data1["bondtypes"] + data2["bondtypes"]
  data["angletypes"] = data1["angletypes"] + data2["angletypes"]
  data["dihedraltypes"] = data1["dihedraltypes"] + data2["dihedraltypes"]

  if box == 1 :
    data["boxx"] = data1["boxx"]
    data["boxy"] = data1["boxy"]
    data["boxz"] = data1["boxz"]
  else :
    data["boxx"] = data2["boxx"]
    data["boxy"] = data2["boxy"]
    data["boxz"] = data2["boxz"]

  data["atype"] = []
  data["atype"].extend(data1["atype"])
  offset = max(data1["atype"])
  for a in data2["atype"] :
    data["atype"].append(a+offset)
  data["xyz"] = np.zeros([data["atoms"],3])
  data["xyz"][:data1["atoms"],:] = data1["xyz"]
  data["xyz"][data1["atoms"]:,:] = data2["xyz"]
  data["q"] = []
  data["q"].extend(data1["q"])
  data["q"].extend(data2["q"])
  data["molecule"] = []
  data["molecule"].extend(data1["molecule"])
  offset = max(data["molecule"])
  for m in data2["molecule"] :
    data["molecule"].append(m+offset)
  if "diameter" in data1 or "diameter" in data2 :
    data["diameter"] = []
    data["density"] = []
    data["mu"] = np.zeros([data["atoms"],3])
    if "diameter" in data1 :
      data["diameter"].extend(data1["diameter"])
      data["density"].extend(data1["density"])
      data["mu"][:data1["atoms"],:] = data1["mu"]
    else :
      data["diameter"].extend([0.0]*data1["atoms"])
      data["density"].extend([1.0]*data1["atoms"])
    if "diamter" in data2 :
      data["diameter"].extend(data2["diameter"])
      data["density"].extend(data2["density"])
      data["mu"][data1["atoms"]:,:] = data2["mu"]
    else :
      data["diameter"].extend([0.0]*data2["atoms"])
      data["density"].extend([1.0]*data2["atoms"])

  if data["bonds"] > 0 :
    data["bonddata"] = np.zeros([data["bonds"],3],dtype=int)
    merge_topol("bond")

  if data["angles"] > 0 :
    data["angledata"] = np.zeros([data["angles"],4],dtype=int)
    merge_topol("angle")

  if data["dihedrals"] > 0 :
    data["dihedraldata"] = np.zeros([data["dihedrals"],5],dtype=int)
    merge_topol("dihedral")

  if data["impropers"] > 0 :
    data["improperdata"] = np.zeros([data["impropers"],5],dtype=int)
    merge_topol("improper")

  return data

def combine_forcefields(forcefield1,forcefield2) :

  forcefield = {}
  forcefield["atomtypes"] = forcefield1["atomtypes"] + forcefield2["atomtypes"]
  forcefield["angletypes"] = forcefield1["angletypes"] + forcefield2["angletypes"]
  forcefield["bondtypes"] = forcefield1["bondtypes"] + forcefield2["bondtypes"]
  forcefield["dihedraltypes"] = forcefield1["dihedraltypes"] + forcefield2["dihedraltypes"]

  forcefield["masses"] = []
  forcefield["masses"].extend(forcefield1["masses"])
  forcefield["masses"].extend(forcefield2["masses"])

  if forcefield["bondtypes"] > 0 :
    forcefield["bondcoeff"] = []
    if forcefield1["bondtypes"] > 0 : forcefield["bondcoeff"].extend(forcefield1["bondcoeff"])
    if forcefield2["bondtypes"] > 0 : forcefield["bondcoeff"].extend(forcefield2["bondcoeff"])

  if forcefield["dihedraltypes"] > 0 :
    forcefield["dihedralcoeff"] = []
    if forcefield1["dihedraltypes"] > 0 : forcefield["dihedralcoeff"].extend(forcefield1["dihedralcoeff"])
    if forcefield2["dihedraltypes"] > 0 : forcefield["dihedralcoeff"].extend(forcefield2["dihedralcoeff"])

  if forcefield["angletypes"] > 0 :
    forcefield["anglecoeff"] = []
    if forcefield1["angletypes"] > 0 : forcefield["anglecoeff"].extend(forcefield1["anglecoeff"])
    if forcefield2["angletypes"] > 0 : forcefield["anglecoeff"].extend(forcefield2["anglecoeff"])
    if "anglepot" in forcefield1 or "anglepot" in forcefield2 :
      forcefield["anglepot"] = []
      if "anglepot" in forcefield1 :
        forcefield["anglepot"].extend(forcefield1["anglepot"])
      else :
        forcefield["anglepot"].extend(["harmonic"]*forcefield1["angletypes"])
      if "anglepot" in forcefield2 :
        forcefield["anglepot"].extend(forcefield2["anglepot"])
      else :
        forcefield["anglepot"].extend(["harmonic"]*forcefield2["angletypes"])

  forcefield["paircoeff"] = np.zeros([forcefield["atomtypes"],forcefield["atomtypes"],2])
  forcefield["pair_an"] = []

  if "paircoeff" in forcefield1 :
    forcefield["paircoeff"][:forcefield1["atomtypes"],:forcefield1["atomtypes"],0] = forcefield1["paircoeff"][:,:,0]
    forcefield["paircoeff"][:forcefield1["atomtypes"],:forcefield1["atomtypes"],1] = forcefield1["paircoeff"][:,:,1]
  elif "ljcoeff" in forcefield1 :
    for i in range(forcefield1["atomtypes"]) :
      for j in range(i,forcefield1["atomtypes"]) :
        eps = np.sqrt(forcefield1["ljcoeff"][i][0]*forcefield1["ljcoeff"][j][0])
        sig = (forcefield1["ljcoeff"][i][1]+forcefield1["ljcoeff"][j][1])/2.0
        forcefield["paircoeff"][i,j,0] = eps
        forcefield["paircoeff"][i,j,1] = sig
  if "pair_an" in forcefield1 :
    for i,a in enumerate(forcefield1["pair_an"]) :
      if len(a) > 2 :
        forcefield["pair_an"].append(a)
      else :
        forcefield["pair_an"].append("From ff 1")
  else :
    for i in range(forcefield1["atomtypes"]) :
      for j in range(i,forcefield1["atomtypes"]) :
        forcefield["pair_an"].append("From ff 1")


  for i in range(forcefield1["atomtypes"])  :
    for j in range(forcefield2["atomtypes"]) :
      forcefield["pair_an"].append("Mixed")

  if "paircoeff" in forcefield2 :
    forcefield["paircoeff"][forcefield1["atomtypes"]:,forcefield1["atomtypes"]:,:] = forcefield2["paircoeff"]
  elif "ljcoeff" in forcefield2 :
    for i in range(forcefield2["atomtypes"]) :
      for j in range(i,forcefield2["atomtypes"]) :
        eps = np.sqrt(forcefield2["ljcoeff"][i][0]*forcefield2["ljcoeff"][j][0])
        sig = (forcefield2["ljcoeff"][i][1]+forcefield2["ljcoeff"][j][1])/2.0
        forcefield["paircoeff"][forcefield1["atomtypes"]+i,forcefield1["atomtypes"]+j,0] = eps
        forcefield["paircoeff"][forcefield1["atomtypes"]+i,forcefield1["atomtypes"]+j,1] = sig
  if "pair_an" in forcefield2 :
    for i,a in enumerate(forcefield2["pair_an"]) :
      if len(a) > 2 :
        forcefield["pair_an"].append(a)
      else :
        forcefield["pair_an"].append("From ff 2")
  else :
    for i in range(forcefield2["atomtypes"]) :
      for j in range(i,forcefield2["atomtypes"]) :
        forcefield["pair_an"].append("From ff 2")

  for i in range(forcefield1["atomtypes"])  :
    for j in range(forcefield2["atomtypes"]) :
      jj = j + forcefield1["atomtypes"]
      eps = np.sqrt(forcefield["paircoeff"][i,i,0]*forcefield["paircoeff"][jj,jj,0])
      sig = (forcefield["paircoeff"][i,i,1]+forcefield["paircoeff"][jj,jj,1])/2.0
      forcefield["paircoeff"][i,forcefield1["atomtypes"]+j,0] = eps
      forcefield["paircoeff"][i,forcefield1["atomtypes"]+j,1] = sig
      forcefield["pair_an"].append("Mixed")

  return forcefield

def read_datafile(filename,data,forcefield,atomtype) :

  data["atoms"] = 0
  data["bonds"] = 0
  data["angles"] = 0
  data["dihedrals"] = 0
  data["impropers"] = 0
  data["atomtypes"] = 0
  data["bondtypes"] = 0
  data["angletypes"] = 0
  data["dihedraltypes"] = 0
  forcefield["atomtypes"] = 0
  forcefield["bondtypes"] = 0
  forcefield["angletypes"] = 0
  forcefield["dihedraltypes"] = 0

  lines = open(filename,'r').readlines()
  data["title"] = lines[0].strip()

  i = 0
  while lines[i].find("zlo") == -1 and i < len(lines)-1:

    i = i + 1

    if lines[i].find(" atoms") > -1 :
      data["atoms"] = int(lines[i].strip().split()[0])
    if lines[i].find(" bonds") > -1 :
      data["bonds"] = int(lines[i].strip().split()[0])
    if lines[i].find(" angles") > -1 :
      data["angles"] = int(lines[i].strip().split()[0])
    if lines[i].find(" dihedrals") > -1 :
      data["dihedrals"] = int(lines[i].strip().split()[0])
    if lines[i].find(" impropers") > -1 :
      data["impropers"] = int(lines[i].strip().split()[0])

    if lines[i].find(" atom types") > -1 :
      data["atomtypes"] = int(lines[i].strip().split()[0])
      forcefield["atomtypes"] = data["atomtypes"]
    if lines[i].find(" bond types") > -1 :
      data["bondtypes"] = int(lines[i].strip().split()[0])
      forcefield["bondtypes"] = data["bondtypes"]
    if lines[i].find(" angle types") > -1 :
      data["angletypes"] = int(lines[i].strip().split()[0])
      forcefield["angletypes"] = data["angletypes"]
    if lines[i].find(" dihedral types") > -1 :
      data["dihedraltypes"] = int(lines[i].strip().split()[0])
      forcefield["dihedraltypes"] = data["dihedraltypes"]

    if lines[i].find("xlo") > -1 :
      cols = lines[i].strip().split()
      data["boxx"] = (float(cols[0]),float(cols[1]))
    if lines[i].find("ylo") > -1 :
      cols = lines[i].strip().split()
      data["boxy"] = (float(cols[0]),float(cols[1]))
    if lines[i].find("zlo") > -1 :
      cols = lines[i].strip().split()
      data["boxz"] = (float(cols[0]),float(cols[1]))

  if data["atoms"] > 0 :
    data["xyz"] = np.zeros([data["atoms"],3])
    data["q"] = np.zeros(data["atoms"])
    data["atype"] = np.zeros(data["atoms"],dtype=int)
    if atomtype != "aa" :
      data["mu"] = np.zeros([data["atoms"],3])
      data["density"] = np.zeros(data["atoms"])
      data["diameter"] = np.zeros(data["atoms"])
      data["molecule"] = np.zeros(data["atoms"],dtype=int)
    else :
      data["molecule"] = np.zeros(data["atoms"],dtype=int)+1

  if data["bonds"] > 0 :
    data["bonddata"] = np.zeros([data["bonds"],3],dtype=int)
  if data["angles"] > 0 :
    data["angledata"] = np.zeros([data["angles"],4],dtype=int)
  if data["dihedrals"] > 0 :
    data["dihedraldata"] = np.zeros([data["dihedrals"],5],dtype=int)
  if data["impropers"] > 0 :
    data["improperdata"] = np.zeros([data["impropers"],5],dtype=int)
  if data["atomtypes"] > 0 :
    forcefield["masses"] = [0]*data["atomtypes"]
    forcefield["ljcoeff"] = [0]*data["atomtypes"]
  if data["bondtypes"] > 0 :
    forcefield["bondcoeff"] = [0]*data["bondtypes"]
  if data["angletypes"] > 0 :
    forcefield["anglecoeff"] = [0]*data["angletypes"]
  if data["dihedraltypes"] > 0 :
    forcefield["dihedralcoeff"] = [0]*data["dihedraltypes"]

  while i < len(lines) :

    if i < len(lines) and lines[i].find("Masses") > - 1 :
      i = i + 2
      for j in range(data["atomtypes"]) :
        forcefield["masses"][j] = float(lines[i].strip().split()[1])
        i = i + 1

    if i < len(lines) and lines[i].find("Pair Coeffs") > - 1 :
      i = i + 2
      for j in range(data["atomtypes"]) :
        forcefield["ljcoeff"][j] = np.array(lines[i].strip().split()[1:],dtype=float)
        i = i + 1

    if i < len(lines) and lines[i].find("Bond Coeffs") > - 1 :
      i = i + 2
      for j in range(data["bondtypes"]) :
        forcefield["bondcoeff"][j] = np.array(lines[i].strip().split()[1:],dtype=float)
        i = i + 1

    if i < len(lines) and lines[i].find("Angle Coeffs") > - 1 :
      i = i + 2
      for j in range(data["angletypes"]) :
        forcefield["anglecoeff"][j] = np.array(lines[i].strip().split()[1:],dtype=float)
        i = i + 1

    if i < len(lines) and lines[i].find("Dihedral Coeffs") > - 1 :
      i = i + 2
      for j in range(data["dihedraltypes"]) :
        forcefield["dihedralcoeff"][j] = np.array(lines[i].strip().split()[1:],dtype=float)
        i = i + 1

    if i < len(lines) and lines[i].find("Bonds") > - 1 :
      i = i + 2
      for j in range(data["bonds"]) :
        data["bonddata"][j,:] = np.array(lines[i].strip().split()[1:],dtype=int)
        i = i + 1

    if i < len(lines) and lines[i].find("Angles") > - 1 :
      i = i + 2
      for j in range(data["angles"]) :
        data["angledata"][j,:] = np.array(lines[i].strip().split()[1:],dtype=int)
        i = i + 1

    if i < len(lines) and lines[i].find("Dihedrals") > - 1 :
      i = i + 2
      for j in range(data["dihedrals"]) :
        data["dihedraldata"][j,:] = np.array(lines[i].strip().split()[1:],dtype=int)
        i = i + 1

    if i < len(lines) and lines[i].find("Atoms") > - 1 :
      i = i + 2
      for j in range(data["atoms"]) :
        cols = lines[i].strip().split()
        if atomtype == "aa" :
          data["xyz"][j,:] = np.array(cols[4:7],dtype=float)
          data["q"][j] = float(cols[3])
          data["atype"][j] = int(cols[2])
        elif atomtype == "cg1" :
          data["xyz"][j,:] = np.array(cols[2:5],dtype=float)
          data["q"][j] = float(cols[6])
          data["atype"][j] = int(cols[1])
          data["mu"][j,:] = np.array(cols[7:10],dtype=float)
          data["diameter"][j] = float(cols[10])
          data["density"][j] = float(cols[11])
          data["molecule"][j] = int(cols[5])
        elif atomtype == "cg2" :
          data["xyz"][j,:] = np.array(cols[2:5],dtype=float)
          data["q"][j] = float(cols[7])
          data["atype"][j] = int(cols[1])
          data["mu"][j,:] = np.array(cols[8:11],dtype=float)
          data["diameter"][j] = float(cols[5])
          data["density"][j] = float(cols[6])
          data["molecule"][j] = int(cols[11])
        i = i + 1

    i = i + 1

def read_forcefield(filename,forcefield) :

  lines = open(filename,'r').readlines()

  forcefield["atomtypes"] = 0
  forcefield["bondtypes"] = 0
  forcefield["angletypes"] = 0
  forcefield["dihedraltypes"] = 0

  i=0
  while i < len(lines) :
    if lines[i].find("mass") == 0 : forcefield["atomtypes"] = forcefield["atomtypes"] + 1
    if lines[i].find("bond_coeff") == 0 : forcefield["bondtypes"] = forcefield["bondtypes"] + 1
    if lines[i].find("angle_coeff") == 0 : forcefield["angletypes"] = forcefield["angletypes"] + 1
    if lines[i].find("dihedral_coeff") == 0 : forcefield["dihedraltypes"] = forcefield["dihedraltypes"] + 1
    i = i + 1

  forcefield["masses"] = np.zeros(forcefield["atomtypes"])
  forcefield["mass_an"] = []
  forcefield["paircoeff"] = np.zeros([forcefield["atomtypes"],forcefield["atomtypes"],2])
  forcefield["pair_style"] = []
  forcefield["pair_an"] = []
  forcefield["bondcoeff"] = np.zeros([forcefield["bondtypes"],2])
  forcefield["bond_an"] = []
  forcefield["anglecoeff"] = np.zeros([forcefield["angletypes"],2])
  if forcefield["angletypes"] > 0 :
    forcefield["angle_an"] = []
    forcefield["anglepot"] = []
  forcefield["dihedralcoeff"] = np.zeros([forcefield["dihedraltypes"],3])
  forcefield["dihedral_an"] = []

  i=0
  while i < len(lines) :

    if i < len(lines) and lines[i].find("mass") == 0 :
      for j in range(forcefield["atomtypes"]) :
        cols = lines[i].strip().split()
        forcefield["masses"][j] = float(cols[2])
        if len(cols) >= 5 :
          forcefield["mass_an"].append(cols[4])
        i = i + 1

    if i < len(lines) and lines[i].find("pair_coeff") == 0 :
      cols = lines[i].strip().split()
      ai = int(cols[1])-1
      aj = int(cols[2])-1
      try :
        test = float(cols[3])
        forcefield["paircoeff"][ai,aj,:] = np.array(cols[3:5],dtype=float)
        if len(cols) >= 7 :
          forcefield["pair_an"].append(cols[6])
        else :
          forcefield["pair_an"].append("")
        try :
          del forcefield["pair_style"]
        except :
          pass
      except :
        forcefield["paircoeff"][ai,aj,:] = np.array(cols[4:6],dtype=float)
        if len(cols) >= 8 :
          forcefield["pair_an"].append(cols[7])
        else :
          forcefield["pair_an"].append("")
        forcefield["pair_style"].append(cols[3])

    if i < len(lines) and lines[i].find("bond_coeff") == 0 :
      for j in range(forcefield["bondtypes"]) :
        forcefield["bondcoeff"][j,:] =  np.array(lines[i].strip().split()[2:4],dtype=float)
        if len(cols) >= 6 :
          forcefield["bond_an"].append(lines[i].strip().split()[5])
        else :
          forcefield["bond_an"].append("")
        i = i + 1

    if i < len(lines) and lines[i].find("angle_coeff") == 0 :
      cols = lines[i].strip().split()
      j = len(forcefield["angle_an"])
      try :
        test = float(cols[2])
        forcefield["anglecoeff"][j,:] =  np.array(cols[2:4],dtype=float)
        del forcefield["anglepot"]
        if len(cols) >= 6 :
          forcefield["angle_an"].append(cols[5])
        else :
          forcefield["angle_an"].append("")
      except :
        forcefield["anglecoeff"][j,:] =  np.array(cols[3:5],dtype=float)
        forcefield["anglepot"].append(cols[2])
        if len(cols) >= 7 :
          forcefield["angle_an"].append(cols[6])
        else :
          forcefield["angle_an"].append("")

    i = i + 1


def write_datafile(filename,data) :

  f = open(filename,"w")

  f.write("%s\n\n"%data["title"])
  f.write("%9d atoms\n"%(data["atoms"]))
  if "bonds" in data :
    f.write("%9d bonds\n"%(data["bonds"]))
  if "angles" in data :
    f.write("%9d angles\n"%(data["angles"]))
  if "dihedrals" in data :
    f.write("%9d dihedrals\n"%(data["dihedrals"]))

  f.write("\n%9d atom types\n"%(data["atomtypes"]))
  if "bondtypes" in data :
    f.write("%9d bond types\n"%(data["bondtypes"]))
  if "angletypes" in data :
    f.write("%9d angle types\n"%(data["angletypes"]))
  if "dihedraltypes" in data :
    f.write("%9d dihedral types\n"%(data["dihedraltypes"]))

  f.write("\n%13.5f%13.5f    xlo xhi\n"%(data["boxx"][0],data["boxx"][1]))
  f.write("%13.5f%13.5f    ylo yhi\n"%(data["boxy"][0],data["boxy"][1]))
  f.write("%13.5f%13.5f    zlo zhi\n"%(data["boxz"][0],data["boxz"][1]))

  f.write("\nAtoms\n\n")
  for i in range(data["atoms"]) :
    if "diameter" in data :
      f.write("%5d%3d%10.4f%10.4f%10.4f%8.3f%8.3f%9.5f%10.4f%10.4f%10.4f%5d\n"%(i+1,data["atype"][i],data["xyz"][i,0],data["xyz"][i,1],data["xyz"][i,2],data["diameter"][i],data["density"][i],data["q"][i],data["mu"][i,0],data["mu"][i,1],data["mu"][i,2],data["molecule"][i]))
    else :
      f.write("%5d%3d%10.4f%10.4f%10.4f%9.5f%5d\n"%(i+1,data["atype"][i],data["xyz"][i,0],data["xyz"][i,1],data["xyz"][i,2],data["q"][i],data["molecule"][i]))

  if "bonds" in data :
    f.write("\nBonds\n\n")
    for i in range(data["bonds"]) :
      f.write("%7d%5d%7d%7d\n"%(i+1,data["bonddata"][i,0],data["bonddata"][i,1],data["bonddata"][i,2]))

  if "angles" in data :
    f.write("\nAngles\n\n")
    for i in range(data["angles"]) :
      f.write("%7d%5d%7d%7d%7d\n"%(i+1,data["angledata"][i,0],data["angledata"][i,1],data["angledata"][i,2],data["angledata"][i,3]))

  if "dihedrals" in data :
    f.write("\nDihedrals\n\n")
    for i in range(data["dihedrals"]) :
      f.write("%7d%5d%7d%7d%7d%7d\n"%(i+1,data["dihedraldata"][i,0],data["dihedraldata"][i,1],data["dihedraldata"][i,2],data["dihedraldata"][i,3],data["dihedraldata"][i,4]))

  f.close()

def write_forcefield(filename,forcefield) :

  f = open(filename,"w")

  f.write("#   atomType  mass \n")
  for i in range(forcefield["atomtypes"]) :
    if "mass_an" in forcefield :
      f.write("mass%5d%9.3f # %s\n"%(i+1,forcefield["masses"][i],forcefield["mass_an"][i]))
    else :
      f.write("mass%5d%9.3f\n"%(i+1,forcefield["masses"][i]))

  f.write("\n# Lennard-Jones coefficients: \n")
  if "pair_style" in forcefield :
    style = "%-25s"%"style"
  else :
    style = ""
  f.write("#           iType jTyp %s eps_ij sig_ij \n"%style)
  k = 0
  for i in range(forcefield["atomtypes"]) :
    for j in range(i,forcefield["atomtypes"]) :
      if "pair_style" in forcefield :
        pair_style = "%-25s"%forcefield["pair_style"][k]
      else :
        pair_style = ""
      if "pair_an" in forcefield and forcefield["pair_an"][k] != "" :
        an = " # %s"%forcefield["pair_an"][k]
      else :
        an = ""
      f.write("pair_coeff%5d%5d %s %9.3f %9.3f %s\n"%(i+1,j+1,pair_style,forcefield["paircoeff"][i,j,0],forcefield["paircoeff"][i,j,1],an))
      k = k + 1

  if "bondtypes" in forcefield :
    f.write("\n# harmonic bond coefficients:\n")
    f.write("#         bondType   K     r0  \n")
    for i in range(forcefield["bondtypes"]) :
      if "bond_an" in forcefield and forcefield["bond_an"][i] != "" :
        f.write("bond_coeff%5d%9.3f%9.3f # %s\n"%(i+1,forcefield["bondcoeff"][i][0],forcefield["bondcoeff"][i][1],forcefield["bond_an"][i]))
      else :
        f.write("bond_coeff%5d%9.3f%9.3f\n"%(i+1,forcefield["bondcoeff"][i][0],forcefield["bondcoeff"][i][1]))

  if "angletypes" in forcefield :
    f.write("\n")
    if "anglepot" in forcefield :
      anglepot = ""
      for i in range(forcefield["angletypes"]) :
        if anglepot != forcefield["anglepot"][i] :
          anglepot = forcefield["anglepot"][i]
          f.write("# %s angle coefficients:\n"%anglepot)
          if anglepot in ["cosine/squared","harmonic"] :
            f.write("#	  angleType	             K    theta0 \n")
          elif anglepot == "dipole" :
            f.write("#	  angleType	   K   gamma0  \n")
      if "angle_an" in forcefield and forcefield["angle_an"][i] != ""  :
        f.write("angle_coeff%5d %15s %9.3f%9.3f # %s\n"%(i+1,anglepot,forcefield["anglecoeff"][i][0],forcefield["anglecoeff"][i][1],forcefield["angle_an"][i]))
      else :
        f.write("angle_coeff%5d %15s %9.3f%9.3f \n"%(i+1,anglepot,forcefield["anglecoeff"][i][0],forcefield["anglecoeff"][i][1]))
    else :
      f.write("# angle coefficients:\n")
      f.write("#	  angleType  K    theta0 \n")
      for i in range(forcefield["angletypes"]) :
        if "angle_an" in forcefield and forcefield["angle_an"][i] != ""  :
          f.write("angle_coeff%5d %9.3f%9.3f # %s\n"%(i+1,forcefield["anglecoeff"][i][0],forcefield["anglecoeff"][i][1],forcefield["angle_an"][i]))
        else :
          f.write("angle_coeff%5d %9.3f%9.3f\n"%(i+1,forcefield["anglecoeff"][i][0],forcefield["anglecoeff"][i][1]))

  if "dihedraltypes" in forcefield :
    f.write("\n# dihedral coefficients:\n")
    for i in range(forcefield["dihedraltypes"]) :
      f.write("dihedral_coeff%5d %9.3f%5d%5d\n"%(i+1,forcefield["dihedralcoeff"][i][0],forcefield["dihedralcoeff"][i][1],forcefield["dihedralcoeff"][i][2]))


  f.close()

if __name__ == '__main__' :

  ff = Includefile(filename=sys.argv[1])
  ff.write(sys.argv[2])
