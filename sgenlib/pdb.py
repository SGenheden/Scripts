# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Classes to read, write and manipulate PDB files

The module contains the following public classes:
  - PDBFile -- the top-level structural class,
    contains chains, residues and atoms
  - Residue -- class to hold a collection of atoms
  - Atom -- class to represent an ATOM or HETATOM record
"""

import sys
import copy

import numpy as np
from scipy.spatial.distance import cdist

codes = {"arg":"R", "his":"H","hid":"H","hie":"H","hip":"H", "lys":"K", "asp":"D", "glu":"E", "ser":"S", "thr":"T", "asn":"N", "gln":"Q", "cys":"C", "gly":"G", "pro":"P", "ala":"A", "val":"V", "ile":"I", "leu":"L", "met":"M", "phe":"F", "tyr":"Y", "trp":"W"}

heavy_aa = {"ALA":["N","CA","CB","C","O"],
"GLY":["N","CA","C","O"],
"SER":["N","CA","CB","OG","C","O"],
"THR":["N","CA","CB","CG2","OG1","C","O"],
"LEU":["N","CA","CB","CG","CD1","CD2","C","O"],
"ILE":["N","CA","CB","CG2","CG1","CD1","C","O"],
"VAL":["N","CA","CB","CG1","CG2","C","O"],
"ASN":["N","CA","CB","CG","OD1","ND2","C","O"],
"GLN":["N","CA","CB","CG","CD","OE1","NE2","C","O"],
"ARG":["N","CA","CB","CG","CD","NE","CZ","NH1","NH2","C","O"],
"HIS":["N","CA","CB","CG","ND1","CE1","NE2","CD2","C","O"],
"HID":["N","CA","CB","CG","ND1","CE1","NE2","CD2","C","O"],
"HIE":["N","CA","CB","CG","ND1","CE1","NE2","CD2","C","O"],
"HIP":["N","CA","CB","CG","ND1","CE1","NE2","CD2","C","O"],
"TRP":["N","CA","CB","CG","CD1","NE1","CE2","CZ2","CH2","CZ3","CE3","CD2","C","O"],
"PHE":["N","CA","CB","CG","CD1","CE1","CZ","CE2","CD2","C","O"],
"TYR":["N","CA","CB","CG","CD1","CE1","CZ","OH","CE2","CD2","C","O"],
"GLU":["N","CA","CB","CG","CD","OE1","OE2","C","O"],
"ASP":["N","CA","CB","CG","OD1","OD2","C","O"],
"LYS":["N","CA","CB","CG","CD","CE","NZ","C","O"],
"PRO":["N","CD","CG","CB","CA","C","O"],
"CYS":["N","CA","CB","SG","C","O"],
"MET":["N","CA","CB","CG","SD","CE","C","O"]}

std_aa_names = {}
std_aa_names["ARG"]="N H CA HA CB HB2 HB3 CG HG2 HG3 CD HD2 HD3 NE HE CZ NH1 HH11 HH12 NH2 HH21 HH22 C O".split()
std_aa_names["HID"]="N H CA HA CB HB2 HB3 CG ND1 HD1 CE1 HE1 NE2 CD2 HD2 C O".split()
std_aa_names["HIE"]="N H CA HA CB HB2 HB3 CG ND1 CE1 HE1 NE2 HE2 CD2 HD2 C O".split()
std_aa_names["HIP"]="N H CA HA CB HB2 HB3 CG ND1 HD1 CE1 HE1 NE2 HE2 CD2 HD2 C O".split()
std_aa_names["HIS"]="N H CA HA CB HB2 HB3 CG ND1 HD1 CE1 HE1 NE2 CD2 HD2 C O".split()
std_aa_names["LYS"]="N H CA HA CB HB2 HB3 CG HG2 HG3 CD HD2 HD3 CE HE2 HE3 NZ HZ1 HZ2 HZ3 C O".split()
std_aa_names["ASP"]="N H CA HA CB HB2 HB3 CG OD1 OD2 C O".split()
std_aa_names["GLU"]="N H CA HA CB HB2 HB3 CG HG2 HG3 CD OE1 OE2 C O".split()
std_aa_names["SER"]="N H CA HA CB HB2 HB3 OG HG C O".split()
std_aa_names["THR"]="N H CA HA CB HB CG2 HG21 HG22 HG23 OG1 HG1 C O".split()
std_aa_names["ASN"]="N H CA HA CB HB2 HB3 CG OD1 ND2 HD21 HD22 C O".split()
std_aa_names["GLN"]="N H CA HA CB HB2 HB3 CG HG2 HG3 CD OE1 NE2 HE21 HE22 C O".split()
std_aa_names["CYS"]="N H CA HA CB HB2 HB3 SG HG C O".split()
std_aa_names["GLY"]="N H CA HA2 HA3 C O".split()
std_aa_names["PRO"]="N CD HD2 HD3 CG HG2 HG3 CB HB2 HB3 CA HA C O".split()
std_aa_names["ALA"]="N H CA HA CB HB1 HB2 HB3 C O".split()
std_aa_names["VAL"]="N H CA HA CB HB CG1 HG11 HG12 HG13 CG2 HG21 HG22 HG23 C O".split()
std_aa_names["ILE"]="N H CA HA CB HB CG2 HG21 HG22 HG23 CG1 HG12 HG13 CD1 HD11 HD12 HD13 C O".split()
std_aa_names["LEU"]="N H CA HA CB HB2 HB3 CG HG CD1 HD11 HD12 HD13 CD2 HD21 HD22 HD23 C O".split()
std_aa_names["MET"]="N H CA HA CB HB2 HB3 CG HG2 HG3 SD CE HE1 HE2 HE3 C O".split()
std_aa_names["PHE"]="N H CA HA CB HB2 HB3 CG CD1 HD1 CE1 HE1 CZ HZ CE2 HE2 CD2 HD2 C O".split()
std_aa_names["TYR"]="N H CA HA CB HB2 HB3 CG CD1 HD1 CE1 HE1 CZ OH HH CE2 HE2 CD2 HD2 C O".split()
std_aa_names["TRP"]="N H CA HA CB HB2 HB3 CG CD1 HD1 NE1 HE1 CE2 CZ2 HZ2 CH2 HH2 CZ3 HZ3 CE3 HE3 CD2 C O".split()
std_aa_names["CYX"]="N H CA HA CB HB2 HB3 SG C O".split()

def make_pdbres(coords,atom_names,res_name,pdbfile) :
  """
  Adds a residue + atoms to a PDBFile structure

  Parameters
  ----------
  coords : Numpy array
    the Cartesian coordinates
  atom_names : list of strings
    the atom names
  res_name : string
    the residue name
  pdbfile : PDBFile object
    the structure to add the residue to
  """
  res = Residue()
  for i,(coord,name) in enumerate(zip(coords,atom_names)) :
    patom = Atom()
    patom.idx = len(pdbfile.atoms)
    patom.hetatm = False
    patom.serial = len(pdbfile.atoms)+1
    patom.name = name
    patom.residue = len(pdbfile.residues)+1
    patom.resname = res_name
    patom.xyz = coord
    patom.x = coord[0]
    patom.y = coord[1]
    patom.z = coord[2]
    res.append(patom)
    pdbfile.atoms.append(patom)
  pdbfile.residues.append(res)

class PDBFile :
  """
  A class to encapsulate atoms, residues and chains from a PDB-file

  Attributes:
  -----------
  atom : list of Atom objects
    the atoms
  residue : list of Residue objects
    the residues
  chains : list of tuples of integers
    the chains
  xyz : numpy array
    the Cartesian coordinates
  filename : string
    the name of the file that was parsed into this object
  charged : list of Residue objects
    the charge residues
  box : numpy array
    the box information
  """
  def __init__(self,filename=None, **kwargs)  :
    self.atoms = []
    self.residues = []
    self.chains = []
    self.xyz = None
    self.filename = filename
    self.charged = None
    self.box = None
    if filename != None :
      self.read(filename, **kwargs)

  def atomindex(self, ambmask) :
      atomnam = ambmask.split('@')[1]
      resname = ambmask.split('@')[0][1:]
      for i, atom in enumerate(self.atoms) :
          if atom.name.strip() == atomnam and atom.resname.strip() == resname :
              return i
      return None

  def extend(self,other) :
    """
    Extend this structure with atoms and residues from another structure

    Parameters
    ----------
    other : PDBFile object
      the structure to extend self with
    """
    for atom in other.atoms :
      self.atoms.append(atom)
    for residue in other.residues :
      self.residues.append(residue)
    self.__parse_chains()

  def extend_residues(self,residues,makecopy=True,dochains=True,resnumber=None) :
    """
    Extend this structure with atoms and residues from a list of residues

    Parameters
    ----------
    residues : list of Residue objects
      the residues to extend self with
    makecopy : boolean, optional
      if to make a copy of the residues or not
    dochains : boolean, optional
      if to make chains
    resnumber : list of integers, optional
      the new residue numbers
    """
    for i,residue in enumerate(residues,1) :
      if makecopy :
        newres = copy.deepcopy(residue)
      else :
        newres = residue
      if resnumber is not None : newres.serial = resnumber[min(i,len(resnumber)-1)]
      self.residues.append(newres)
      for atom in newres.atoms :
        self.atoms.append(atom)
        if resnumber is not None : atom.residue = resnumber[min(i,len(resnumber)-1)]
      if dochains : self.__parse_chains()

  def renumber(self,doatoms=True,doresidues=True) :
    """
    Renumber atoms and residues from 1

    Parameters
    ----------
    doatoms : boolean, optional
      if to renumber atoms
    doresidues : boolean, optional
      if to renumber residues
    """
    if doresidues :
      for i,residue in enumerate(self.residues,1) :
        residue.serial = i
        for atom in residue.atoms : atom.residue = i
    if doatoms :
      for i,atom in enumerate(self.atoms,1) :
        atom.serial = i

  def reorder(self,resnames,atomnames=[]) :
    """
    Reorder the residues and atoms

    Parameters
    ----------
    resnames : list of string
      the order of the residues
    atomnames : list of list of strings
      the atom names of each sorted residue
    """
    new_residues = []
    for nam in resnames :
      i = 0
      while i < len(self.residues) :
        if self.residues[i].resname.strip() == nam.strip() :
          new_residues.append(self.residues.pop(i))
        i += 1
    self.residues = new_residues

    if len(atomnames) != len(self.residues) :
      atomnames = [None]*len(self.residues)

    self.atoms = []
    for res,nams in zip(self.residues,atomnames) :
      if nams is not None : res.reorder(nams)
      self.atoms.extend(res.atoms)

    self.renumber()

  def split_residue(self,index,npieces) :
    """
    Split a residue into several individual residues
    """

    new_len = len(self.residues[index].atoms)/npieces
    for j in range(1,npieces) :
      new_res = Residue()
      for i in range(new_len) :
        new_res.append(self.residues[index].atoms.pop(new_len))
      self.residues.append(new_res)

  def cterminal(self) :
    """
    Return a list of the indices of residues at the end of chains
    """
    return [chain[1] for chain in self.chains]

  def nterminal(self) :
    """
    Return a list of the indices of residues at the start of chains
    """
    return [chain[0] for chain in self.chains]

  def charged_residues(self) :
    """
    Produce, store and return a list of charged residues (including C- and N-terminals)
    """
    if not self.residues : return []
    if self.charged != None : return self.charged

    self.charged = []
    taken = [False]*len(self.residues)
    for i,residue in enumerate(self.residues) :
      if residue.resname in ["ASP","GLU","LYS","ARG","HIS","HIP"] :
        self.charged.append(residue)
        taken[i] = True
    # Add N and C-terminal residues that not has been taken already
    for chain in self.chains :
      if not taken[chain[0]] : self.charged.append(self.residues[chain[0]])
      if not taken[chain[1]] : self.charged.append(self.residues[chain[1]])
    return self.charged

  def update_xyz(self,xyz) :
     """
     Update the Cartesian coordinates of all atoms
     """
     self.xyz = np.array(xyz,copy=True)
     for atom,coord in zip(self.atoms,xyz) :
       atom.x = coord[0]
       atom.y = coord[1]
       atom.z = coord[2]
       atom.xyz = np.array(coord,copy=True)

  def read(self,filename,gro=False, **kwargs) :
    """
    Read a PDB-file or GRO-file and parse into atoms, residues and chains

    Parameters
    ----------
    filename : string
      the name of the file to parse
    gro : boolean, optional
      if to force parsing of GRO structure
    """
    self.filename = filename
    if filename[-3:].lower() == "gro" or gro :
      self.__parse_gro_records()
    else :
      self.__parse_records()
    self.__parse_residues(**kwargs)
    self.__parse_chains()

  def write(self,filename=None,ter=False,add_extra=None) :
    """
    Write the structure to a file

    Parameters
    ----------
    filename : string, optional
      the name of the file to write to
    ter : boolean, optional
      if to add TER records
    add_extra : function, optional
      function to add extra, non ATOM, HETATM and TER records
    """
    if filename is None :
      filename = self.filename
    if filename[-3:].lower() == "gro" :
      self.write_gro(filename)
    else:
      with open(filename,"w") as f :
        self.write_to(f,ter=ter,add_extra=add_extra)

  def write_to(self,f,ter=False,add_extra=None) :
    """
    Write the structure to a file object

    Parameters
    ----------
    f : File object
      the file object to write to
    ter : boolean, optional
      if to add TER records
    add_extra : function, optional
      function to add extra, non ATOM, HETATM and TER records
    """
    if self.box is not None :
      f.write("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2fP 1\n"%(self.box[0],self.box[1],self.box[2],90.0,90.0,90.0))
    # If TER record not wanted, the just print out each residue
    if not ter :
      for residue in self.residues :
        f.write(residue.__str__())
    # If TER record is wanted, then we need to do something more
    else :
      # Print out each chain in turn, followed by a TER
      for chain in self.chains :
        for i in range(chain[0],chain[1]+1) :
          f.write(self.residues[i].__str__())
        f.write("TER\n")
      # Then write out each HET residues, followed by a TER
      if len(self.chains) > 0 :
        for i in range(self.chains[-1][1]+1,len(self.residues)) :
          if self.residues[i].hidden : continue
          f.write(self.residues[i].__str__())
          f.write("TER\n")
      else :
        for residue in self.residues :
          f.write(residue.__str__())
          f.write("TER\n")
    if add_extra : add_extra(self,f)

  def write_gro(self,filename=None) :
    """
    Write the structure to a GRO-file

    filename : string, optional
      the name of the file to write to
    """
    if filename == None :
      filename = self.filename
    f = open(filename,"w")
    f.write("Written by pdb.py\n")
    n = 0
    for atom in self.atoms :
      if not atom.hidden : n += 1
    f.write("%8d\n"%n)
    aserial = 1
    for atom in self.atoms :
      if atom.hidden : continue
      resid = (atom.residue if atom.residue <= 99999 else atom.residue - 99999)
      serial = (aserial if aserial <= 99999 else aserial - 99999)
      f.write("%5d%5s%5s%5d%8.3f%8.3f%8.3f\n"%(resid,atom.resname,atom.name,serial,atom.x/10.0,atom.y/10.0,atom.z/10.0))
      aserial += 1
    if self.box is not None :
      f.write("%8.3f%8.3f%8.3f\n"%(self.box[0]/10.0,self.box[1]/10.0,self.box[2]/10.0))
    else :
       f.write("%8.3f%8.3f%8.3f\n"%(0.0,0.0,0.0))

  def __parse_chains(self) :
    """
    Find start and end of chains from an array of Residue objects
    This only identify chains of ATOM records, not HETATM
    """
    self.chains = []
    first = 0
    last = -1
    prev = self.residues[0].chain
    for i,residue in enumerate(self.residues) :
      if i > 0 and residue.atoms[0].hetatm :
        last = i-1
        self.chains.append((first,last))
        break
      if residue.chain != prev :
        last = i-1
        self.chains.append((first,last))
        first = i
        last = -1
      prev = residue.chain

    if not self.residues[-1].atoms[0].hetatm :
      last = i-1
      self.chains.append((first,last))

  def __parse_residues(self, renumber=True) :
    """
    Produce an array of Residue objectes from an array of Atom objects
    """
    residues = []
    self.residues.append(Residue(idx=0,atom=self.atoms[0]))

    for atom in self.atoms[1:] :
      if atom.residue != self.residues[-1].serial or atom.resname != self.residues[-1].resname  or atom.chain != self.residues[-1].chain :
        self.residues.append(Residue(idx=len(self.residues),atom=atom))
      else :
        self.residues[-1].append(atom)

    if renumber :
      for i,residue in enumerate(self.residues,1):
        residue.set_serial(i)

  def __parse_records(self) :
    """
    Read ATOM/HETATM records from PDB file and create an array of Atom objects
    and a numpy array of coordinates
    """
    self.atoms = []
    self.xyz = []
    with open(self.filename,"r") as f :
      line = f.readline()
      while line :
        if line[0:6] in ["ATOM  ","HETATM"] :
          atom = Atom(record=line,idx=len(self.atoms))
          self.atoms.append(atom)
          self.xyz.append([atom.x,atom.y,atom.z])
        elif line[0:6] == "CRYST1" :
          self.box = np.array(line.strip().split()[1:4],dtype=float)
        line = f.readline()
    self.xyz = np.array(self.xyz)

  def __parse_gro_records(self) :
    """
    Read ATOM/HETATM records from GRO file and create an array of Atom objects
    and a numpy array of coordinates
    """
    self.atoms = []
    self.xyz = []
    lines = open(self.filename,"r").readlines()
    for line in lines[2:-1] :
      atom = Atom(idx=len(self.atoms))
      atom.readGRO(line)
      self.atoms.append(atom)
      self.xyz.append([atom.x,atom.y,atom.z])

    self.xyz = np.array(self.xyz)
    self.box = np.array(lines[-1].strip().split(),dtype=float)*10.0

class Atom :
  """
  A class to encapsulate an atom from a PDB-file

  Attributes
  ----------
  idx : integer
    the serial number in the PDB-file atom list
  hetatm : boolean
    if this a HETATM record
  serial : integer
    the serial number as read from the PDB-file
  name : string
    the atom name
  altloc : string
    the alternative location
  resname : string
    the residue name
  chain : string
    the chain identifier
  residue : integer
    the residue serial number
  insertion : string
    the insertion code
  x, y, z : float
    the Cartesian coordinates
  occupancy : float
    the occupancy
  bfactor : float
    the bfactor
  term : string
    the rest of the record
  xyz : numpy array
    the Cartesian coordinates
  hidden : boolean
    indicates if this should be written to file
  pqratom : boolean
    indicates if this should be written in a free-format
  """
  def __init__(self,idx=None,record=None) :
    self.idx = idx
    self.hetatm = False
    self.serial = 0
    self.name = ""
    self.altloc = ""
    self.resname = ""
    self.chain = ""
    self.residue = 0
    self.insertion = ""
    self.x = 0.0
    self.y = 0.0
    self.z = 0.0
    self.occupancy = 0.0
    self.bfactor = 0.0
    self.term = ""
    self.xyz = None
    self.hidden = False
    self.pqratom = False
    if record :
      self.readRecord(record)

  def __str__(self) :
    """
    Produces a string representation of the atom, i.e. an ATOM/HETATM PDB record
    """
    if self.hidden : return ""
    record = {False:"ATOM  ",True:"HETATM"}
    name = self.name
    if len(self.name) > 4 :
      name = self.name[:4]
    resname  = self.resname
    if len(self.resname) > 3 :
      resname = self.resname[:3]
    if self.term == "" or self.term[-1] != "\n" : self.term = self.term+ "\n"
    if self.residue > 9999 :
      resnr = self.residue - 9999
    else :
      resnr = self.residue
    if self.serial > 99999 :
      serial = self.serial - 99999
    else :
      serial = self.serial
    if not self.pqratom :
        return "%6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f%s"%(record[self.hetatm],serial,name,self.altloc,resname,self.chain,resnr,self.insertion,self.x,self.y,self.z,self.occupancy,self.bfactor,self.term)
    else :
        return "%6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%8.4f%8.4f%s"%("ATOM  ",serial,name,self.altloc,resname,self.chain,resnr,self.insertion,self.x,self.y,self.z,self.occupancy,self.bfactor,self.term)

  def readGRO(self,record) :
    """
    Read a record from a GRO-file
    """
    self.serial = int(record[15:20].strip())
    self.name = record[10:15].strip()
    #if len(self.name) > 4 :
    #  self.name = self.name[:4]
    #elif len(self.name) < 4 :
    #  self.name = "%4s"%self.name
    self.resname = record[5:10].strip()
    #if len(self.resname) > 3 :
    #  self.resname = self.resname[:3]
    #elif len(self.resname) < 3 :
    #  self.resname = "%3s"%self.resname
    self.residue = int(record[:5].strip())
    self.x = float(record[20:28].strip())*10
    self.y = float(record[28:36].strip())*10
    self.z = float(record[36:44].strip())*10
    self.xyz = np.array([self.x,self.y,self.z])
    self.hetatm = not self.resname.upper() in heavy_aa.keys()
    self.term = "\n"
    self.altloc = ""
    self.chain = ""
    self.insertions = ""
    self.occupancy = 0.0
    self.bfactor = 0.0
    self.hidden = False

  def readRecord(self,record) :
    """
    Read a PDB HETATM or ATOM record
    """
    self.hetatm = record[0:6] == "HETATM"
    test = record[6:11].strip()
    if test.find("*") > -1 :
      self.serial = -1
    else :
      self.serial = int(record[6:11].strip())
    self.name = record[12:16]
    self.altloc = record[16]
    self.resname = record[17:20]
    self.chain = record[21]
    self.residue = int(record[22:27].strip())
    self.insertion = record[26]
    self.x = float(record[30:38].strip())
    self.y = float(record[38:46].strip())
    self.z = float(record[46:54].strip())
    try :
      self.occupancy = float(record[54:60].strip())
    except :
      self.occupancy = 0.0
    try :
      self.bfactor = float(record[60:66].strip())
    except :
      self.bfactor = 0.0
    try :
      self.term = record[66:]
    except :
      self.term = "\n"
    self.xyz = np.array([self.x,self.y,self.z])
    self.hidden = False

  def distance2(self,atom) :
    """
    Calculates the squared distances to another atom
    """
    return (self.x-atom.x)**2+(self.y-atom.y)**2+(self.z-atom.z)**2

  def mass(self) :
    """
    Returns the mass of this atom
    """
    masses = {"h":1.00800,"c":12.01100,"n":14.00700,"o":15.99940,"p":30.97400,"s":32.066}
    if self.element() not in masses :
      raise Exception("Do not know the mass of element %s"%self.element())
    return masses[self.element()]

  def element(self) :
    """
    Returns the element of this atom
    ONLY works for single character elements
    """
    name = self.name.strip().lower()
    try :
      trial = int(name[0])
      name = name[1:]
    except :
      pass
    return name[0]

  def set_xyz(self,xyz) :
    """
    Sets the Cartesian coordinates
    """
    self.x = xyz[0]
    self.y = xyz[1]
    self.z = xyz[2]
    self.xyz = np.array(xyz)

class Residue :
  """
  Class to encapsulate a collection of Atom objects

  Attributes
  ----------
  atoms : list of Atom objects
    the atoms
  idx : integer
    the serial number of this residue in the PDBFile structure
  hidden : boolean
    if this residue should be written to file
  serial : integer
    the serial number as read from disc
  resname : string
    the name of this residue
  chain : string
    the chain identifier
  """
  def __init__(self,idx=None,atom=None) :
    self.atoms = []
    self.idx = idx
    self.hidden = False
    if atom :
      self.atoms.append(atom)
      self.serial = atom.residue
      self.resname = atom.resname
      self.chain = atom.chain
      self.hidden = atom.hidden
    else :
      self.serial = 1
      self.resname = ""
      self.chain   = ""

  def __str__(self) :
    """
    Produces a string representation of all the atoms in the residues
    """
    if self.hidden : return ""
    return "".join([atom.__str__() for atom in self.atoms])

  def append(self,atom) :
    """
    Append another Atom object to the collection
    """
    self.atoms.append(atom)
    if len(self.atoms) == 1 :
      self.serial = atom.residue
      self.resname = atom.resname
      self.chain = atom.chain
      self.hidden = atom.hidden
  def atom_by_name(self, name) :
      """
      Return an atom in this residue by its name
      """
      for atom in self.atoms:
          if atom.name.strip() == name :
              return atom
      return None
  def index_by_name(self, name) :
       """
       Return an atom index in this residue by its name
       """
       for i, atom in enumerate(self.atoms):
           if atom.name.strip() == name :
               return i
       return None
  def check_heavy(self) :
    """
    Check if an amino acid residue is complete
    """
    # Return if this is an empty collection or this residue is not an amino acid
    if len(self.atoms) == 0 or self.resname.upper() not in heavy_aa.keys() : return None,None

    found = {}
    for aname in heavy_aa[self.resname.upper()] :
      found[aname] = 0
    found["OXT"] = 0
    extra = []
    for atom in self.atoms :
      aname = atom.name.strip().upper()
      if aname != "H" :
        if aname in heavy_aa[self.resname.upper()] or aname == "OXT" :
          found[aname] = found[aname] + 1
        else :
          extra.append(aname)
    notfound = []
    for aname in heavy_aa[self.resname.upper()] :
      if found[aname] == 0 : notfound.append(aname)
    return notfound,extra

  def set_hidden(self,hidden) :
    """
    Updates the hidden flag
    """
    self.hidden = hidden
    for atom in self.atoms :
      atom.hidden = hidden

  def set_resname(self,resname) :
    """
    Updates the residue name
    """
    #resname = resname.upper()
    self.resname = resname
    for atom in self.atoms :
      atom.resname = resname

  def set_serial(self,serial) :
    """
    Updates the serial number
    """
    self.serial = serial
    for atom in self.atoms :
      atom.residue = serial

  def update_xyz(self,xyz) :
    """
    Update the Cartesian coordinates of this residue
    """
    for i,atom in enumerate(self.atoms) :
      atom.x = xyz[i,0]
      atom.y = xyz[i,1]
      atom.z = xyz[i,2]
      atom.xyz = np.array(xyz[i,:],copy=True)

  def distance2(self,residue,xyz) :
    """
    Calculates the minimum squared distance to another residue

    Attributes
    ----------
    residue : Residue object
      the other residue object
    xyz : numpy array
      all Cartesian coordinates of the PDBFile
    """
    xyz1 = xyz[self.atoms[0].idx:self.atoms[-1].idx+1,:]
    xyz2 = xyz[residue.atoms[0].idx:residue.atoms[-1].idx+1,:]
    mindist = cdist(xyz1,xyz2,"sqeuclidean").min()
    return mindist

  def within(self,residue,xyz,cutoff) :
    """
    Return whether this residues is within a cutoff of another residue
    """
    return self.distance2(residue,xyz) <= cutoff*cutoff

  def collect(self,operation) :
    """
    Perform an operation on the collection

    Can do center of mass, masses and xyz
    """
    vector = np.zeros([len(self.atoms),3])
    operation = operation.lower()
    for i,atom in enumerate(self.atoms) :
      if operation in ["centerofmass","xyz"] :
        vector[i,:] = atom.xyz
      elif operation == "masses" :
        vector[i,0] = atom.mass()

    if operation == "centerofmass" :
      return np.mean(vector,axis=0)
    elif operation == "masses" :
      return vector[:,0]
    return vector

  def reorder(self, atomnames) :
    """
    Reorder the atoms according to a given list of atom names
    """
    new_atoms = []
    for nam in atomnames :
      i = 0
      while i < len(self.atoms) :
        if self.atoms[i].name.strip() == nam.strip() :
          new_atoms.append(self.atoms.pop(i))
        i += 1
    self.atoms = new_atoms

  def rename(self,atomnames) :
    """
    Rename atoms according to the dictionary
    """
    for atom in self.atoms :
      aname = atom.name.strip()
      if not aname in atomnames :
        raise Exception("Could not find %s in dictionary"%aname)
      else :
        if atomnames[aname] == "*" :
          atom.hidden = True
        else :
          atom.name = atomnames[aname]
