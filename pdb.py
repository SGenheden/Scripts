import sys
import numpy as np
from scipy.spatial.distance import cdist

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

#
# A class to encapsulate atoms, residues and chains from a PDB-file
#
class PDBFile :
  def __init__(self,filename=None)  :
    self.atoms = []
    self.residues = []
    self.chains = []
    self.xyz = None
    self.filename = filename
    self.charged = None
    self.box = None
    if filename != None :
      self.read(filename)
  
  def extend(self,other) :
    for atom in other.atoms :
      self.atoms.append(atom)
    for residue in other.residues :
      self.residues.append(residue)
    self.__parse_chains()    
  
  # Return a list of the indices of residues at the end of chains    
  def cterminal(self) :
    return [chain[1] for chain in self.chains]
    
  # Return a list of the indices of residues at the start of chains
  def nterminal(self) :
    return [chain[0] for chain in self.chains]
    
  # Produce, store and return a list of charged residues (including C- and N-terminals)    
  def charged_residues(self) :
  
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
    
  #
  # Read a PDB-file or GRO-file and parse into atoms, residues and chains
  #
  def read(self,filename,gro=False) :
    self.filename = filename
    if not gro :
      self.__parse_records()
    else :
      self.__parse_gro_records()
    self.__parse_residues()
    self.__parse_chains()

  #
  # Write all residues to file
  #
  def write(self,filename=None,ter=False,add_extra=None) :
    if filename == None :
      filename = self.filename
    f = open(filename,"w")
    if self.box != None :
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
          f.write(self.residues[i].__str__())
          f.write("TER\n")
      else :
        for residue in self.residues :
          f.write(residue.__str__())
          f.write("TER\n")
    if add_extra : add_extra(self,f)
    f.close()

  #
  #
  #
  def write_gro(self,filename=None) :
    if filename == None :
      filename = self.filename
    f = open(filename,"w")
    f.write("Written by pdb.py\n")
    f.write("%8d\n"%len(self.atoms))
    for atom in self.atoms :
      f.write("%5d%5s%5s%5d%8.3f%8.3f%8.3f\n"%(atom.residue,atom.resname,atom.name,atom.serial,atom.x/10.0,atom.y/10.0,atom.z/10.0))
    if self.box != None :
      f.write("%8.3f%8.3f%8.3f\n"%(self.box[0]/10.0,self.box[1]/10.0,self.box[2]/10.0))
    else :
       f.write("%8.3f%8.3f%8.3f\n"%(0.0,0.0,0.0))

#          1         2         3         4
#0123456789012345678901234567890123456789012345678901234567890123456789
# 5248SOL     H233024   4.444   0.859   0.806  0.6749 -0.4063  0.6694
    
  #
  # Find start and end of chains from an array of Residue objects 
  # This only identify chains of ATOM records, not HETATM 
  #
  def __parse_chains(self) :
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
  #
  # Produce an array of Residue objectes from an array of Atom objects
  #
  def __parse_residues(self) :

    residues = []
    self.residues.append(Residue(idx=0,atom=self.atoms[0]))
    
    for atom in self.atoms[1:] :
      if atom.residue != self.residues[-1].serial or atom.resname != self.residues[-1].resname  or atom.chain != self.residues[-1].chain :
        self.residues.append(Residue(idx=len(self.residues),atom=atom))
      else :
        self.residues[-1].append(atom)
  #
  # Read ATOM/HETATM records from PDB file and create an array of Atom objects
  # and a numpy array of coordinates
  #
  def __parse_records(self) :

    self.atoms = []
    self.xyz = []
    f = open(self.filename,"r")
    line = f.readline()
    while line :
      if line[0:6] in ["ATOM  ","HETATM"] :
        atom = Atom(record=line,idx=len(self.atoms))
        self.atoms.append(atom)
        self.xyz.append([atom.x,atom.y,atom.z])
      elif line[0:6] == "CRYST1" :
        self.box = np.array(line.strip().split()[1:4],dtype=float)
      line = f.readline()
    f.close()

    self.xyz = np.array(self.xyz)

  #
  # Read ATOM/HETATM records from GRO file and create an array of Atom objects
  # and a numpy array of coordinates
  #
  def __parse_gro_records(self) :

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
#
# A class to encapsulate an atom from a PDB-file
#
class Atom :
  def __init__(self,idx=None,record=None,copy=None) :
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
    if record :
      self.readRecord(record)
    if copy :
      self.copy(copy)
  #
  # Produces a string representation of the atom, i.e. an ATOM/HETATM PDB record
  #
  def __str__(self) :
    if self.hidden : return ""
    record = {False:"ATOM  ",True:"HETATM"}
    name = self.name
    if len(self.name) > 4 :
      name = self.name[:4]
    resname  = self.resname
    if len(self.resname) > 3 :
      resname = self.resname[:3]
    if self.term == "" or self.term[-1] != "\n" : self.term = self.term+ "\n"
    return "%6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f%s"%(record[self.hetatm],self.serial,name,self.altloc,resname,self.chain,self.residue,self.insertion,self.x,self.y,self.z,self.occupancy,self.bfactor,self.term)    
  #
  # Read a record from a gro-file
  #
#          1         2         3         4
#01234567890123456789012345678901234567890123456789
#    1CHOL   ROH    1  -1.050   1.030   2.808
  def readGRO(self,record) :
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
  #
  # Read a PDB HETATM or ATOM record
  #
  def readRecord(self,record) :
    self.hetatm = record[0:6] == "HETATM"
    self.serial = int(record[6:11].strip())
    self.name = record[12:16]
    self.altloc = record[16]
    self.resname = record[17:20]
    self.chain = record[21]
    self.residue = int(record[22:26].strip())
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
  # 
  # Copy all properties of the atom
  #
  def copy(self,copy) :
    self.idx = copy.idx
    self.hetatm = copy.hetatm
    self.serial = copy.serial
    self.name = copy.name
    self.altloc = copy.altloc
    self.resname = copy.resname
    self.chain = copy.chain
    self.residue = copy.residue
    self.insertion = copy.insertion
    self.x = copy.x
    self.y = copy.y
    self.z = copy.z
    self.occupancy = copy.occupancy
    self.bfactor = copy.bfactor
    self.term = copy.term
    self.xyz = np.array(copy.xyz,copy=True)
    self.hidden = copy.hidden
  #
  # Calculates the squared distances to another atom
  #
  def distance2(self,atom) :
    return (self.x-atom.x)**2+(self.y-atom.y)**2+(self.z-atom.z)**2
  #
  #
  #
  def mass(self) : 
    masses = {"h":1.00800,"c":12.01100,"n":14.00700,"o":15.99940,"p":30.97400,"s":32.066}
    if self.element() not in masses :
      raise Exception("Do not know the mass of element %s"%self.element())
    return masses[self.element()]
  #
  # 
  def element(self) :
    name = self.name.strip().lower()
    try :
      trial = int(name[0])
      name = name[1:]
    except :
      pass
    return name[0]

#
# Class to encapsulate a list of Atom objects
#
class Residue :
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
      self.serial = ""
      self.resname = ""
      self.chain   = ""
  #
  # Produces a string representation of all the atoms in the residues
  #
  def __str__(self) :
    if self.hidden : return ""
    return "".join([atom.__str__() for atom in self.atoms])
  #
  # Append another Atom object to the collection
  #
  def append(self,atom) :
    self.atoms.append(atom)
    if len(self.atoms) == 1 :
      #self.atoms.append(atom)
      self.serial = atom.residue
      self.resname = atom.resname
      self.chain = atom.chain
      self.hidden = atom.hidden  
  #
  # Check if an amino acid residue is complete
  #
  def check_heavy(self) :
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
  #
  # Make a copy of the residue
  #
  def copy(self,residue) :
    self.atoms = []
    for atom in residue.atoms :
      self.append(Atom(copy=atom))
  # 
  # Updates the residue name
  #  
  def set_resname(self,resname) :
    resname = resname.upper()
    if len(resname) > 3 : 
      resname = resname[:3]
    elif len(resname) < 3 :
      resname = "%3s"%resname
      
    self.resname = resname
    for atom in self.atoms :
      atom.resname = resname

  def set_serial(self,serial) :
    self.serial = serial
    for atom in self.atoms :
      atom.residue = serial

  #
  # Calculates the minimum squared distance to another residue 
  #
  def distance2(self,residue,xyz) :
    xyz1 = xyz[self.atoms[0].idx:self.atoms[-1].idx+1,:]
    xyz2 = xyz[residue.atoms[0].idx:residue.atoms[-1].idx+1,:]
    mindist = cdist(xyz1,xyz2,"sqeuclidean").min()     
    return mindist
  #
  # Return whether this residues is within a cutoff of another residue
  #
  def within(self,residue,xyz,cutoff) :
    return self.distance2(residue,xyz) <= cutoff*cutoff
  #
  #
  #
  def collect(self,operation) :
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
  #
  #
  #
  def update_xyz(self,xyz) :
    for i,atom in enumerate(self.atoms) :
      atom.x = xyz[i,0]
      atom.y = xyz[i,1]
      atom.z = xyz[i,2]
      atom.xyz = np.array(xyz[i,:],copy=True)
