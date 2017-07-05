# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Classes to read and manipulate some Gromacs files
"""

import os

from . import geo


#
# Class to hold an AtomType record
#
class AtomType :
  def __init__(self,record=None) :
    self.name = ""
    self.atnum = 0
    self.mass = 0.0
    self.charge = 0.0
    self.sigma = 0.0
    self.epsilon = 0.0
    self.record = ""
    if record != None : self.read(record)
  def read(self,record) :
    self.record = record
    cols = record.strip().split()
    # for Martini, no general
    try :
      dummy = float(cols[3])
    except :
      cols.insert(1,0)
    self.name = cols[0]
    self.atnum = cols[1]
    try :
      self.mass = float(cols[2])
    except :
      self.mass = 0.0
    if len(cols) > 3 :
      self.charge = float(cols[3])
      offset = 1
      try :
        self.sigma = float(cols[5])
        offset = 0
      except :
        offset = 1
      self.sigma = float(cols[5+offset])
      self.epsilon = float(cols[6+offset])
  def __str__(self) :
    return "%s %s %.2f %.2f %.2f %.2f"%(self.name,self.atnum,self.mass,self.charge,self.sigma,self.epsilon)

#
# Class to hold an PairType record
#
class PairType :
  def __init__(self,record=None) :
    self.atomj = ""
    self.atomi = ""
    self.sigma = 0.0
    self.epsilon = 0.0
    self.record = ""
    if record != None : self.read(record)
  def read(self,record) :
    self.record = record
    cols = record.strip().split()
    if not cols[2].isdigit() or int(cols[2]) != 1 : return
    self.atomi = cols[0]
    self.atomj = cols[1]
    self.sigma = float(cols[3])
    self.epsilon = float(cols[4])
  def __str__(self) :
    return "%s %s %.2E %.2E"%(self.atomj,self.atomi,self.sigma,self.epsilon)

#
# Class to hold a Connectivity (bond,angle,dihedral) type record
#
class ConnectivityType :
  def __init__(self,record=None,natoms=None) :
    self.record = ""
    self.atoms = []
    self.func = 0
    self.params = []
    self.addparams = []
    self.addfunc = []
    if record != None : self.read(record,natoms)
  def read(self,record,natoms) :
    self.record = record
    i = record.find(";")
    if i > -1 : record = record[:i]
    cols = record.strip().split()
    self.atoms = cols[:natoms]
    try :
      self.func = int(cols[natoms])
    except :
      pass
    self.params = []
    for i in cols[natoms+1:] :
      try :
        ii  = float(i)
        self.params.append(ii)
      except :
        self.params.append(i)
  def __str__(self) :
    atomstr = " ".join(self.atoms)
    paramstr = " ".join(["%.2f"%i for i in self.params])
    if len(self.addparams) > 0 :
      for params in self.addparams :
        paramstr = paramstr + "\n\t %s"%" ".join(["%.2f"%i for i in params])
    return "%s %d %s"%(atomstr,self.func,paramstr)

#
# Class to hold an atom record
#
class MolAtom :
  def __init__(self,record=None) :
    self.id = 0
    self.type = ""
    self.residue = 0
    self.resname = ""
    self.name = ""
    self.cgrp = 0
    self.charge = 0.0
    self.mass = 0.0
    if record != None : self.read(record)
  def read(self,record) :
    cols = record.strip().split()
    self.id = int(cols[0])
    self.type = cols[1]
    self.residue = int(cols[2])
    self.resname = cols[3]
    self.name = cols[4]
    self.cgrp = int(cols[5])
    self.charge = float(cols[6])
    if len(cols) > 7 :
      try :
        self.mass = float(cols[7])
      except :
        pass
  def __str__(self) :
    return "%6d %8s %6d %6s %6s %4d %10.4f %10.6f"%(self.id,self.type,self.residue,self.resname,self.name,self.cgrp,self.charge,self.mass)

#
# Class to hold a molecular topology (bond,angle,dihedral)
#
class Connectivity :
  def __init__(self,record=None,natoms=None) :
    self.atoms = []
    self.func = 0
    self.params = []
    self.addparams = []
    self.addfunc = []
    if record != None : self.read(record,natoms)
  def read(self,record,natoms) :
    i = record.find(";")
    if i > -1 : record = record[:i]
    cols = record.strip().split()
    self.atoms = [int(i) for i in cols[:natoms]]
    try :
        self.func = int(cols[natoms])
    except :
        pass
    try :
      self.params = [float(i) for i in cols[natoms+1:]]
    except :
      self.params = []
  def __str__(self) :
    atomstr = " ".join(["%d "%i for i in self.atoms])
    paramstr = " ".join(["%.2f"%i for i in self.params])
    return "%s %d %s"%(atomstr,self.func,paramstr)

#
# Class to read a molecular type
#
class MoleculeType :
  def __init__(self) :
    self.name = None
    self.atoms = []
    self.bonds = []
    self.angles = []
    self.dihedrals = []
    self.pairs = []
  # Read several atom records
  def __read_atoms(self,lines,lineidx) :
    while lineidx < len(lines) and lines[lineidx].find("[") == -1 :
      if not (lines[lineidx][0] in [";","#"] or len(lines[lineidx].strip()) < 2) :
        self.atoms.append(MolAtom(record=lines[lineidx]))
      lineidx = lineidx + 1
    return lineidx
  # Read several topology records
  def __read_connectivity(self,lines,lineidx,container,natoms) :
    while lineidx < len(lines) and lines[lineidx].find("[") == -1 :
      if not (lines[lineidx][0] in [";","#"] or len(lines[lineidx].strip()) < 2) :
        container.append(Connectivity(record=lines[lineidx],natoms=natoms))
      lineidx = lineidx + 1
    return lineidx
  # Parse a full molecular type from a file
  def read(self,lines,lineidx) :

    while lines[lineidx][0] == ";" or len(lines[lineidx])  < 2 : lineidx = lineidx + 1
    self.name = lines[lineidx].strip().split()[0]
    
    while lineidx < len(lines) and lines[lineidx].find("[ moleculetype ]") == -1 :
      if lines[lineidx].find("[ atoms ]") > -1 or lines[lineidx].find("[atoms]") > -1:
        lineidx = self.__read_atoms(lines,lineidx+1)
      elif lines[lineidx].find("[ bonds ]") > -1 or lines[lineidx].find("[bonds]") > -1:
        lineidx = self.__read_connectivity(lines,lineidx+1,self.bonds,2)
      elif lines[lineidx].find("[ pairs ]") > -1 or lines[lineidx].find("[pairs]") > -1:
        lineidx = self.__read_connectivity(lines,lineidx+1,self.pairs,2)
      elif lines[lineidx].find("[ angles ]") > -1 or lines[lineidx].find("[angles]") > -1:
        lineidx = self.__read_connectivity(lines,lineidx+1,self.angles,3)
      elif lines[lineidx].find("[ dihedrals ]") > -1 or lines[lineidx].find("[dihedrals]") > -1:
        lineidx = self.__read_connectivity(lines,lineidx+1,self.dihedrals,4)
      else :
        lineidx = lineidx + 1
    return lineidx
  def __str__(self) :
    str = self.name+"\n"
    for atom in self.atoms :
      str = str + "\t%s\n"%atom.__str__()
    str = str + "\nBonds:\n"
    for bond in self.bonds :
      str = str + "\t%s\n"%bond.__str__()
    str = str + "\nAngles:\n"
    for angle in self.angles :
      str = str + "\t%s\n"%angle.__str__()
    str = str + "\nDihedrals:\n"
    for dihedral in self.dihedrals :
      str = str + "\t%s\n"%dihedral.__str__()
    return str

#
# Super class to hold everything parsed from a top-file
#
class TopFile :
  def __init__(self,filename=None) :
    self.atomtypes = []
    self.pairtypes = []
    self.bondtypes = []
    self.angletypes = []
    self.dihedraltypes = []
    self.moleculetypes = []
    self.name = ""
    self.molecules = {}
    self.mollist = []
    self.specialLJ = 0.0
    self.specialCoul = 0.0
    if filename != None : self.read(filename)
  def __include(self,filename) :
    fullpath = os.path.dirname(os.path.abspath(filename))
    lines = open(filename,"r").readlines()
    lineidx = 0
    while lineidx < len(lines) :
      if lines[lineidx].find("#include") == 0:
        self.__include(os.path.join(fullpath,lines[lineidx].strip().split()[1].replace('"',"")))
        lineidx = lineidx + 1
      elif lines[lineidx].find("[ system ]") == 0 :
        lineidx = lineidx + 1
        while lines[lineidx][0] == ";" or len(lines[lineidx])  < 2 : lineidx = lineidx + 1
        self.name = lines[lineidx].strip().split()[0]
      elif lines[lineidx].find("[ defaults ]") == 0 :
        lineidx = lineidx + 1
        while lines[lineidx][0] == ";" or len(lines[lineidx])  < 2 : lineidx = lineidx + 1
        cols = lines[lineidx].strip().split()
        if len(cols) > 2 :
          self.specialLJ = float(cols[3])
          self.specialCoul = float(cols[4])
        else :
          self.specialLJ = 1.0
          self.specialCoul = 1.0
      elif lines[lineidx].find("[ moleculetype ]") == 0 :
        self.moleculetypes.append(MoleculeType())
        lineidx = self.moleculetypes[-1].read(lines,lineidx+1)
      elif lines[lineidx].find("[ atomtypes ]") == 0 :
        lineidx = self.__read_records(lines,lineidx+1,AtomType,self.atomtypes)
      elif lines[lineidx].find("[ pairtypes ]") == 0 or lines[lineidx].find("[ nonbond_params ]") == 0:
        lineidx = self.__read_records(lines,lineidx+1,PairType,self.pairtypes)
      elif lines[lineidx].find("[ bondtypes ]") == 0 :
        lineidx = self.__read_records(lines,lineidx+1,ConnectivityType,self.bondtypes,extra=2)
      elif lines[lineidx].find("[ angletypes ]") == 0 :
        lineidx = self.__read_records(lines,lineidx+1,ConnectivityType,self.angletypes,extra=3)
      elif lines[lineidx].find("[ dihedraltypes ]") == 0 :
        lineidx = self.__read_records(lines,lineidx+1,ConnectivityType,self.dihedraltypes,extra=4)
      elif lines[lineidx].find("[ molecules ]") == 0 :
         lineidx = self.__read_molecules(lines,lineidx+1)
      else :
        lineidx = lineidx + 1
  def __read_molecules(self,lines,lineidx) :
    while lineidx < len(lines) and lines[lineidx].find("[ ") == -1 :
      if not (lines[lineidx][0] == ";" or len(lines[lineidx]) < 2) :
        cols = lines[lineidx].strip().split()
        self.molecules[cols[0]] = int(cols[1])
        self.mollist.append(cols[0])
      lineidx = lineidx + 1
    return lineidx
  def __read_records(self,lines,lineidx,classtype,container,extra=None) :
    while lineidx < len(lines) and lines[lineidx].find("[ ") == -1 :
      if not (lines[lineidx][0] == ";" or len(lines[lineidx].strip()) < 2) :
        if extra == None :
          container.append(classtype(lines[lineidx]))
        else :
          container.append(classtype(lines[lineidx],extra))
      lineidx = lineidx + 1
    return lineidx
  def read(self,filename) :
    self.__include(filename)
    #self.__pair_dihedrals()
  def __pair_dihedrals(self) :
    taken = [False]*len(self.dihedraltypes)
    for i,dihedral1 in enumerate(self.dihedraltypes) :
      if taken[i] or dihedral1.func != 9 : continue
      atomstr1 = " ".join(dihedral1.atoms)
      for j,dihedral2 in enumerate(self.dihedraltypes[i+1:]) :
        if taken[i+j+1] or dihedral1.func != dihedral2.func : continue
        atomstr2 =  " ".join(dihedral2.atoms)
        if atomstr1 == atomstr2 :
          dihedral1.addparams.append(dihedral2.params)
          dihedral1.addfunc.append(dihedral2.func)
          taken[i+j+1] = True
    new_dihedraltypes = []
    for i,dihedral in enumerate(self.dihedraltypes) :
      if not taken[i] : new_dihedraltypes.append(dihedral)
    self.dihedraltypes = new_dihedraltypes
  def reduce2RB(self) :
    new_dihedraltypes = []
    for i,dihedral in enumerate(self.dihedraltypes) :
      if dihedral.func != 9 :
        new_dihedraltypes.append(dihedral)
      else :
        forces = [dihedral.params[1]]
        multiplicities = [int(dihedral.params[2])]
        phases = [dihedral.params[0]]
        for params in dihedral.addparams :
          forces.append(params[1])
          multiplicities.append(int(params[2]))
          phases.append(params[0])
        new_dihedraltypes.append(dihedral)
        new_dihedraltypes[-1].addparams = []
        new_dihedraltypes[-1].func = 3
        new_dihedraltypes[-1].params = geo.multi2RB(forces,multiplicities,phases)
    self.dihedraltypes = new_dihedraltypes

    for mol in self.moleculetypes :
      for dihedral in mol.dihedral :
        if dihedral.func == 9 : dihedral.func = 3

#
# Classes to pair conformation (PBD objects) with topology (Gmx objects)
#

class ConfTopAtom :
  def __init__(self,top,conf) :
    self.top = top
    self.conf = conf

class ConfTopResidue :
  def __init__(self,topmol,atom=None) :
    self.atoms = []
    self.topmol = topmol
    if atom != None : atoms.append(atom)

class ConfTopMol :
  def __init__(self,top,conf) :
    self.atoms = []
    self.residues = []
    self.match(top,conf)
  def match(self,top,conf) :
    for residue in conf.residues :
      # First find a matching residue
      resname = residue.resname.strip().split()
      topmol = None
      for mol in top.moleculetypes :
        for atom in mol.atoms :
          if atom.resname.strip().split() == resname :
            topmol = mol
            break
        if topmol != None : break
      self.residues.append(ConfTopResidue(topmol))
      # Now match atoms in the residue
      for atom in residue.atoms :
        atomname = atom.name.strip().split()
        for atom2 in topmol.atoms :
          if atom2.name.strip().split() == atomname :
            self.atoms.append(ConfTopAtom(atom2,atom))
            self.residues[-1].atoms.append(self.atoms[-1])
