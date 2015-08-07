# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to prepare a PDB-file semi-automatically

The program will download the PDB file and metadata automatically
from the PDB database. It will describe the most important
features about the structure. It will then prompt the user
for removal of chains, alternative conformations, hetero residues
and the protonation state of histidines.

Examples
--------
  pdb_prep.py 18IL

"""

import sys
import copy
import urllib
import os.path
import xml.etree.ElementTree as et
import ConfigParser
import logging

import numpy as np
from scipy.spatial.distance import cdist

from sgenlib.pdb import *


# ------------------
# Utility routines
# ------------------

def _setup_logger(filename) :
  """
  Setup logging system
  
  Uses the standard logging module with two handlers,
  one that prints everything above DEBUG to standard output
  and one that prints everything, even DEBUG to a file
  
  If filename is None no file handler is created
  
  Parameters
  ----------
  filename : string, optional
    the filename of the string
    
  Returns
  -------
  a reference to the created logger
  """
  logger = logging.getLogger('pdb_prep')
  logger.setLevel(logging.DEBUG)
  formatter1 = logging.Formatter('%(message)s')
  console = logging.StreamHandler()
  console.setLevel(logging.INFO)
  console.setFormatter(formatter1)
  logger.addHandler(console)
  if filename is not None :
    formatter2 = logging.Formatter('%(levelname)s : %(message)s')
    logfile = logging.FileHandler(filename,mode="w")
    logfile.setLevel(logging.DEBUG)
    logfile.setFormatter(formatter2)
    logger.addHandler(logfile)
  return logger

def _new_section(wait=True) :
  """
  Write out a delimiter and wait for user response
  
  Parameters
  ----------
  wait : boolean
    if to wait for the user
  """
  if wait : raw_input("\nPress enter to continue")
  logger.info("#"*80)
  logger.info("")

def _fetch(url,localname,label) :
  """
  Download a file from the Internet
  
  Parameters
  ----------
  url : string
    the URL to download
  localname : string
    the location to download it to
  label : string
    a label to indicate for the user

  Raises
  ------
  Exception
    the URL could not be retrieved
  """
  try : 
    logger.info("Fetching %s..."%label)
    path, header = urllib.urlretrieve(url,localname)
  except :
    msg = "Could not retrieve %s!"%label
    logger.error(msg)
    raise Exception(msg)

# ----------------------------------------------
# Routines to read and write selections to disc
# ----------------------------------------------

def _read_sel(sel,pid) :
  """
  Read an ini-style file from disc if it exists otherwise set default values
  
  Parameters
  ----------
  sel : dictionary
    the selections
  pid : string
    the PDB identifier
  """
  parser = ConfigParser.SafeConfigParser()
  parser.read("%s_sel.dat"%pid)
  
  # Set None as default
  sel["pdb"]  = None
  sel["meta"] = None
  sel["ignore_chains"] = None
  sel["ignore_atoms"] = None
  sel["cyscmd"] = None
  sel["his"] = None
  sel["ignore_het"] = None
  sel["buried"] = None
  
  if not parser.has_section("selections") : return 
  
  for option in ["pdb","meta","ignore_chains","ignore_atoms","ignore_het","cyscmd","his"] :
    if parser.has_option("selections",option) :
      if option in ["ignore_chains","ignore_het","his"] :
        sel[option] = parser.get("selections",option).split()
      else :
        sel[option] = parser.get("selections",option)
          
  if sel["ignore_atoms"] is not None : sel["ignore_atoms"] = map(int,sel["ignore_atoms"].split()) 

def _write_sel(sel,pid) :
  """
  Write an ini-style file to disc with user selections
  
  Parameters
  ----------
  sel : dictionary
    the selections
  pid : string
    the PDB identifier 
  """
  parser = ConfigParser.SafeConfigParser()

  parser.add_section('selections')
  # String selections
  for option in ["pdb","meta","cyscmd"] :
    if sel[option] is not None : parser.set("selections",option,sel[option])
  # List of string selections
  for option in ["ignore_chains","ignore_het","his"] :
    if sel[option] is not None : parser.set("selections",option," ".join(sel[option]))
  # List of int selection
  if sel["ignore_atoms"] is not None : parser.set("selections","ignore_atoms"," ".join("%d"%a for a in sel["ignore_atoms"]))
  
  with open("%s_sel.dat"%pid,"w") as f :
    parser.write(f)



# ------------------------------------
# Main routines to prepare a PDB-file
# -------------------------------------


def _do_describe(pdbfile,pid,filename) :
  """
  Print metadata and statistics for a PDB-file
  
  Parameters
  ----------
  meta : dictionary 
    the meta data parsed from a PDB XML-file
  pdbfile : PDBFile object
    the pdb structure
  """
  
  try :  
    xml = et.parse(filename)
  except :
    msg = "Could not parse file retrieved from the PDB!"
    logger.error(msg)
    raise Exception(msg)
  root = xml.getroot()
  roottag = root.tag
  ns = roottag[:-len("datablock")]
  
  meta = {} 
  meta["id"] = pid
  try :
    meta["resolution"] = root.findall("%sreflnsCategory/%sreflns[1]/%sd_resolution_high"%(ns,ns,ns) )[0].text
  except :
    meta["resolution"] = "0"
  meta["date0"] = root.findall("%sdatabase_PDB_revCategory/%sdatabase_PDB_rev[1]/%sdate"%(ns,ns,ns) )[0].text
  meta["date"] = root.findall("%sdatabase_PDB_revCategory/%sdatabase_PDB_rev[last()]/%sdate"%(ns,ns,ns) )[0].text
  meta["title"] = root.findall("%sstructCategory/%sstruct[1]/%stitle"%(ns,ns,ns) )[0].text
  try :
    meta["rfree"] = root.findall("%srefineCategory/%srefine[1]/%sls_R_factor_R_free"%(ns,ns,ns) )[0].text
  except :
    meta["rfree"] = "0"
  meta["molecules"] = []
  for ent in root.findall("%sentityCategory/%sentity"%(ns,ns) ) : 
    meta["molecules"].append(ent.findall("%spdbx_description"%ns)[0].text)
  try :
    meta["src_common"] = root.findall("%sentity_src_genCategory/%sentity_src_gen[1]/%sgene_src_common_name"%(ns,ns,ns) )[0].text
  except :
    meta["src_common"] = None
  try :
    meta["src"] = root.findall("%sentity_src_genCategory/%sentity_src_gen[1]/%spdbx_gene_src_scientific_name"%(ns,ns,ns) )[0].text
  except :
    meta["src"] = ""
  try :
    meta["ph"] = root.findall("%sexptl_crystal_growCategory/%sexptl_crystal_grow[1]/%spH"%(ns,ns,ns) )[0].text
  except :
    meta["ph"] = "0.0"
  
  logger.info("Meta data for PDB file %s"%meta["id"].upper())
  logger.info("Title: %s"%meta["title"])
  logger.info("Date of release: %s\tLast modified: %s"%(meta["date0"],meta["date"]))
  logger.info("Resolution: %s\tRfree: %s\tpH: %s"%(meta["resolution"],meta["rfree"],meta["ph"]))
  if meta["src_common"] :
    logger.info("Source: %s (%s) "%(meta["src_common"],meta["src"]))
  else :
    logger.info("Source: %s "%meta["src"])
  logger.info("Molecules: %s"%", ".join(meta["molecules"]))
  logger.info("")

  stdaa = ["ARG","HIS","LYS","ASP","GLU","SER","THR","ASN","GLN","CYS","GLY","PRO","ALA","VAL","ILE","LEU","MET","PHE","TYR","TRP"]
  
  logger.info("Atoms: %5d\tResidues: %5d"%(len(pdbfile.atoms),len(pdbfile.residues)))
  nwat = naa = nhet = narg = nlys = nhis = nasp = nglu = noxt = 0
  for residue in pdbfile.residues :
    if residue.resname in stdaa : naa = naa + 1
    if residue.resname == "HOH" : nwat = nwat + 1
    if residue.atoms[0].hetatm : nhet = nhet + 1
    if residue.resname == "ARG" :  narg = narg + 1
    if residue.resname == "LYS" : nlys = nlys + 1
    if residue.resname == "HIS" : nhis = nhis + 1
    if residue.resname == "ASP" : nasp = nasp + 1
    if residue.resname == "GLU" : nglu = nglu + 1
  for atom in pdbfile.atoms :
    if atom.name.strip().upper() == "OXT" : noxt = noxt + 1
  logger.info("Standard amino-acids: %4d\tWater molecules: %4d\tHetero residues: %4d"%(naa,nwat,nhet-nwat))
  logger.info("ARG: %4d\tLYS: %4d"%(narg,nlys)) 
  logger.info("ASP: %4d\tGLU: %4d"%(nasp,nglu))
  logger.info("HIS: %4d\tOXT: %4d"%(nhis,noxt))
  logger.info("Charged residues sum to %de"%(narg+nlys-nasp-nglu))
  
  logger.info("")
  allok = True
  for residue in pdbfile.residues :
    if residue.resname.upper() in heavy_aa.keys() :
      notfound,extra= residue.check_heavy()
      if len(notfound) > 0 or len(extra) > 0 : allok = False
      if len(notfound) > 0 :
        logger.info("Missing atom(s) for %s%d: (%s)"%(residue.resname.capitalize(),residue.serial,", ".join(notfound)))
      if len(extra) > 0 :
        logger.info("Extra atom(s) for %s%d: (%s)"%(residue.resname.capitalize(),residue.serial,", ".join(extra)))
  if allok : logger.info("All heavy amino-acid atoms found")
  
  logger.info("")
  logger.info("%5s %8s %8s"%("Chain","Start   ","End     "))
  for chain in pdbfile.chains :
    res0 = pdbfile.residues[chain[0]]
    res1 = pdbfile.residues[chain[1]]
    logger.info("%5s %3s%-5d %3s%-5d"%(res0.chain,res0.resname.capitalize(),res0.serial,res1.resname.capitalize(),res1.serial))

  logger.info("")
  allok = True
  for chain in pdbfile.chains :
    for i in range(chain[0],chain[1]) :
      if pdbfile.residues[i].serial != pdbfile.residues[i+1].serial -1 :
        res0 = pdbfile.residues[i]
        res1 = pdbfile.residues[i+1]
        logger.info("There is a gap in chain %s between %3s%d and %3s%d"%(res0.chain,res0.resname.capitalize(),res0.serial,res1.resname.capitalize(),res1.serial))
        allok = False
  if allok : logger.info("All chains are intact")
  
  return filename

def _do_buried(pdbfile,pid,selection) :
  """
  Find buried charges and ionic pairs
  """
  
  def find_residues() :
    
    buried = []
    threshold = 15.5*15.5
    neighcut = 575
    for residue in pdbfile.charged_residues() :
      nneigh = 0
      xyz0 = np.zeros(3)
      n0 = 0
      for atom in residue.atoms :
        if (residue.resname == "ASP" and atom.name.strip().upper() in ["OD1","OD2"]) or \
           (residue.resname == "GLU" and atom.name.strip().upper() in ["OE1","OE2"]) or \
           (residue.resname == "LYS" and atom.name.strip().upper() == "NZ") or \
           (residue.resname == "ARG" and atom.name.strip().upper() == "CZ") or \
           (residue.resname in ["HIS","HID","HIE"] and atom.name.strip().upper() in ["CG","ND1","CD2","NE2","CE1"]) or \
           (residue.idx in pdbfile.nterminal() and atom.name.strip().upper() == "N") or \
           (residue.idx in pdbfile.cterminal() and atom.name.strip().upper() in ["O","OXT"]) :
           xyz0 = xyz0 + atom.xyz
           n0 = n0 + 1
      if n0 > 0 :
        xyz0 = xyz0 / float(n0)
        for residue2 in pdbfile.residues :
          if residue2.idx != residue.idx and residue2.resname != "HOH" and residue.within(residue2,pdbfile.xyz,threshold) :
             for atom in residue2.atoms :
               if atom.name.strip().lower()[0] != "H" and ((atom.xyz-xyz0)**2).sum() < threshold : nneigh = nneigh + 1
        if nneigh >= neighcut :
          logger.info("Buried charge: %s%d - %.2f"%(residue.resname.capitalize(),residue.serial,float(nneigh)/float(neighcut)*100.0))
          buried.append(residue)
    return buried
  
  buried = find_residues()
  
  logger.info("")
  pairs = {}
  for bresidue in buried :
    for residue in pdbfile.charged_residues() :
      if bresidue.idx != residue.idx : 
        for a1 in bresidue.atoms :
          name1 = a1.name.strip().upper()
          if (not name1[0] == "C") and (not name1 in ["O","N"]) :
            for a2 in residue.atoms :
              name2 = a2.name.strip().upper() 
              if (not name2[0] == "C") and (not name2 in ["O","N"]) and a1.distance2(a2) < (3.5*3.5) :
                 if bresidue.idx < residue.idx :
                   pairs["Possible ionic pair: %s%d - %s%d"%(bresidue.resname.capitalize(),bresidue.serial,residue.resname.capitalize(),residue.serial)] = 1
                 else :
                   pairs["Possible ionic pair: %s%d - %s%d"%(residue.resname.capitalize(),residue.serial,bresidue.resname.capitalize(),bresidue.serial)] = 1
  for pair in pairs.keys() :
    logger.info(pair)
  return None

def _do_cys(pdbfile,pid,user_sel) :
  """
  Find CYS-CYS bridges and ask the user if they should be written to a leap command-file
  """
  cys = []
  for residue in pdbfile.residues :
    if residue.resname.upper()[:2] == "CY" :
      cys.append(residue)

  logger.info("Identified CYS-CYS bridges:")
  threshold = 2.50*2.50
  bridges = []
  for i,cys1 in enumerate(cys) :
    for atom1 in cys1.atoms :
      if atom1.name.strip() == "SG" :
        for cys2 in cys[i+1:] :
          for atom2 in cys2.atoms :
            if atom2.name.strip() == "SG" :
              d = atom1.distance2(atom2)
              if d <= threshold :
                logger.info("%3s%-5d - %3s%-5d\t%.3f A"%(cys1.resname.capitalize(),cys1.serial,cys2.resname.capitalize(),cys2.serial,np.sqrt(d)))
                bridges.append((cys1,cys2))
  logger.info("Found %d CYS residues and %s CYS-CYS bridges"%(len(cys),len(bridges)))

  cyscmd = ""
  if len(bridges) > 0 :
    if user_sel == None :
      sel = raw_input("Do you want to write the bridges to a leap command-file (yes, no)? ")
    else :
      if len(user_sel) > 0 :
        sel = "yes"
      else :
        sel = "no"
    if sel.lower() == "yes" or (len(sel)==1 and sel[0].lower() == "y"): 
      cyscmd = "%s_cyx.dat"%pid.lower() 
      with open(cyscmd,"w") as f :
        for b in bridges :
          f.write("bond x.%d.SG x.%d.SG\n"%(b[0].serial,b[1].serial))
      logger.info("Wrote leap commands to %s"%cyscmd)
      
    for bridge in bridges :
      bridge[0].set_resname("CYX")
      bridge[1].set_resname("CYX")
    logger.info("")
    logger.info("Changed all CYS involved in bridges to CYX")
  return cyscmd

def _do_chains(pdbfile,pid,user_sel) :
  """
  Prompt the user to remove chains
  """
  ignore = []
  chainid = [pdbfile.residues[chain[0]].chain.upper() for chain in pdbfile.chains]
  logger.info("The PDB-file contains %d chains (%s). Do you want to remove some chain(s)?"%(len(pdbfile.chains),",".join(chainid)))
  if user_sel == None :
    sel = raw_input("Remove the following chains: ")
  else :
    sel = " ".join(user_sel)
  if sel == "" : 
    logger.info("Keeping all chains")
  else :
    for s in sel.strip().split() :
      if s.upper() in chainid : 
        chidx = chainid.index(s)
        for residue in pdbfile.residues[pdbfile.chains[chidx][0]:pdbfile.chains[chidx][1]+1] :
          residue.set_hidden(True)
        ignore.append(s.upper())
    logger.info("Ignoring the chain(s) %s"%",".join(ignore))
  return ignore

def _do_conformations(pdbfile,pdbid,user_sel) :
  """
  Prompt the user for removal of alternative conformations
  """
  has_low = []
  logger.info("Residues with occupancy numbers less than 1.0:")
  for i,residue in enumerate(pdbfile.residues) :
    for atom in residue.atoms :
      if atom.occupancy < 1.0 :
        logger.info("%s%d\t%.2f"%(residue.resname.capitalize(),residue.serial,atom.occupancy))
        has_low.append(i)
        break
        
  logger.info("")
  logger.info("Residues with alternative locations:")
  for residue in pdbfile.residues :
    for atom in residue.atoms :
      if atom.altloc != " " :
        logger.info("%s%d"%(residue.resname.capitalize(),residue.serial))
        break

  if user_sel == None :
    sel = raw_input("\nDo you want to remove conformations with lower occupancy (yes,no)? ")
  else :
    if len(user_sel) > 0 :
      sel = "yes"
    else :
      sel = "no"
  ignore = []
  if sel.lower() == "yes" or (len(sel)==1 and sel[0].lower() == "y"):  
    for i in has_low :
      minocc = 1.0
      minloc = "B"
      for atom in pdbfile.residues[i].atoms :
        minocc = min(minocc,atom.occupancy) 
        if minocc == atom.occupancy : minloc = atom.altloc
      if minocc == 0.5 : # Remove B conformation
        logger.info("Removing B conformation for %s%d"%(pdbfile.residues[i].resname.capitalize(),pdbfile.residues[i].serial))
        for atom in pdbfile.residues[i].atoms :
          if atom.altloc == "B" : 
            ignore.append(atom.idx)
            atom.hidden = True
      else : 
        logger.info("Removing %s conformation for %s%d"%(minloc,pdbfile.residues[i].resname.capitalize(),pdbfile.residues[i].serial))
        for atom in pdbfile.residues[i].atoms :
          if atom.occupancy == minocc : 
            ignore.append(atom.idx)
            atom.hidden = True
    logger.info("Removing %d atoms in total"%len(ignore))
  else :
    if len(has_low) > 0 :
      logger.info("Keeping all conformations. Note that this might cause problems.")
    else :
      logger.info("Keeping all conformations.")
  return ignore
  
def _do_het(pdbfile,pid,user_sel) :
  """
  Find hetero residues (excluding water) and ask the user if some of them should be removed
  """
  
  het = []
  hetres = []
  ignore = []
  for residue in pdbfile.residues :
    if residue.atoms[0].hetatm and residue.resname.upper() != "HOH" :
      het.append(residue)
      hetres.append(residue.resname.strip().capitalize())
      
  logger.info("The PDB-file contains the following hetero residues (%s)."%",".join(hetres))
  logger.info("Do you want to remove some of them?")
  if user_sel == None :
    sel = raw_input("Remove the following residues: ")
  else :
    sel = " ".join(user_sel)
  if sel == "" : 
    logger.info("Keeping all residues")
  else :
    for s in sel.strip().split() :
      if s.capitalize() in hetres : ignore.append(s.capitalize())
    logger.info("Ignoring the residue(s) %s"%",".join(ignore))
    for residue in het :
      if residue.resname.strip().capitalize() in ignore : residue.set_hidden(True)
  return ignore

def _do_his(pdbfile,pid,user_sel) :
  """
  Prompt the user for histidine protonation state
  """

  def add_his_hydrogen(n,c1,c2,name) :
    """
    Create histidine protons
    """
    h = copy.deepcopy(n)
    d1 = 2.00*n.xyz-(c1.xyz+c2.xyz)
    ld1 = 1.0300
    f = ld1 / np.sqrt(np.sum(d1**2))
    d1 = d1*f
    h.xyz = n.xyz+d1
    h.x   = n.x+d1[0]
    h.y   = n.y+d1[1]
    h.z   = n.z+d1[2]
    h.name = name
    return h
    
  def angle(v1,v2) :
    """
    Calculates the angles between two vectors
    """
    l1 = np.sqrt(np.sum(v1**2))
    l2 = np.sqrt(np.sum(v2**2))
    a = np.sum(np.multiply(v1,v2))/(l1*l2)
    if a > 1 :
      a  = 0.0
    elif a < -1 :
      a = 180.0
    else :
      a = np.arccos(a)*(180.0/np.pi)
    return a

  def find_hbond(h,n,c1,c2,residue,presidue,donors,acceptors,anyinter,waters) :
    """
    Find h-bond partners to histidine protons
    """
    d1 = h.xyz - n.xyz
    for i,atom in enumerate(residue.atoms) :
      d3 = atom.xyz - h.xyz
      dish = np.sqrt(np.sum(d3**2))
      disn = np.sqrt(n.distance2(atom))
      if atom.name.strip()[0] != "C" :
        d3 = -d3
        ang1 = angle(d1,d3)
      
        if len(residue.atoms) > 2 :
          if residue.atoms[i-1].name.strip()[0] == "C" :
            k = i-1
          else :
            k = i-2
          if k < 0 :
            d4 = residue.atoms[-abs(k)].xyz-atom.xyz
          else :
            d4 = residue.atoms[k].xyz-atom.xyz
          ld4 = np.sqrt(np.sum(d4**2))
          ang2 = angle(d3,d4)
        else :
          ang2 = 180.0

        if dish < 3.0 and ang1 > 90 and ang2 > 90 :
          aname = atom.name.strip().upper()
          rname = atom.resname.upper()
          if aname == "N" :
            donors.append([atom,dish,disn,ang1,ang2])
          elif aname == "O" and rname != "HOH" :
            acceptors.append([atom,dish,disn,ang1,ang2])
          elif rname in ["TRP","ARG","LYS"] :
            donors.append([atom,dish,disn,ang1,ang2])
          elif rname in ["MET","ASP","GLU"] :
            acceptors.append([atom,dish,disn,ang1,ang2])
          elif rname in ["GLN","ASN"] :
            if aname[0] == "O" :
              acceptors.append([atom,dish,disn,ang1,ang2])
            elif aname[0] == "N" :
              donors.append([atom,dish,disn,ang1,ang2])
          elif rname in ["THR","SER","TYR"] :
            anyinter.append([atom,dish,disn,ang1,ang2])          
        elif dish < 3.0 :
          anyinter.append([atom,dish,disn,ang1,ang2])
        elif dish < 2.5 and atom.resname.upper() == "HOH" :
          waters.append([atom,dish,disn,ang1,0.0])   
  
  his = []
  for residue in pdbfile.residues :
    if residue.resname.upper()[:2] == "HI" :
      his.append(residue)  

  threshold = 8.0
  threshold2 = threshold*threshold
  for i,h in enumerate(his) :
    logger.info("")
    logger.info("-- His%d --"%h.serial)
    f = open("his%d.pdb"%(i+1),"w")
    f.write(h.__str__())
    f.write("TER\n")
    cg  = None 
    nd1 = None
    ce1 = None
    cd2 = None
    ne2 = None
    for atom in h.atoms :
      if atom.name.strip().upper() == "CG" :  cg  = atom
      if atom.name.strip().upper() == "ND1" : nd1 = atom
      if atom.name.strip().upper() == "CE1" : ce1 = atom
      if atom.name.strip().upper() == "CD2" : cd2 = atom
      if atom.name.strip().upper() == "NE2" : ne2 = atom
    hd1 = add_his_hydrogen(nd1,cg,ce1,"HD1")
    he2 = add_his_hydrogen(ne2,cd2,ce1,"HE2")
    
    nneigh = 0
    sumch = 0.0
    hid_hbonds = {"Donor: ":[],"Acceptor: ":[],"Any: ":[],"Water: ":[]}
    hie_hbonds = {"Donor: ":[],"Acceptor: ":[],"Any: ":[],"Water: ":[]}
    presidue = None
    chres = []
    for residue in pdbfile.residues :
      if residue.resname != h.resname and residue.serial != h.serial and h.within(residue,pdbfile.xyz,threshold) :
        f.write(residue.__str__())
        f.write("TER\n")
        for atom in residue.atoms :
          if ce1.distance2(atom) < threshold2 :
            nneigh = nneigh + 1
            if (residue.resname.upper() == "ASP" and atom.name.upper().strip() == "CG") or (residue.resname.upper() == "GLU" and atom.name.upper().strip() == "CD") : 
              chres.append("%s%d"%(residue.resname.capitalize(),residue.serial))
              sumch = sumch -1
            if (residue.resname.upper() == "ARG" and atom.name.upper().strip() == "CZ") or (residue.resname.upper() == "LYS" and atom.name.upper().strip() == "NZ") :
              chres.append("%s%d"%(residue.resname.capitalize(),residue.serial))
              sumch = sumch +1
        find_hbond(hd1,nd1,cg,ce1,residue,presidue,hid_hbonds["Donor: "],hid_hbonds["Acceptor: "],hid_hbonds["Any: "],hid_hbonds["Water: "])
        find_hbond(he2,ne2,cd2,ce1,residue,presidue,hie_hbonds["Donor: "],hie_hbonds["Acceptor: "],hie_hbonds["Any: "],hie_hbonds["Water: "])
      presidue = residue    
    f.close()
    logger.info(" Atoms within 8 A of CE1: %3d\tSolvent exposure of CE1 : %d%%\t"%(nneigh,int((1-float(nneigh-15)/115.0)*100.0)))
    logger.info("Charge within 8 A of CE1: %3d\t (%s)"%(sumch,", ".join(chres)))
    logger.info("")
    if len(hid_hbonds["Donor: "]) > 0 or len(hid_hbonds["Acceptor: "]) > 0 or len(hid_hbonds["Any: "]) > 0 or len(hid_hbonds["Water: "]) > 0 :
      logger.info("**HD1**    \tResidue   Atom\t %8s %8s %8s %8s"%("r(X-H)","r(X-N)","v(NHX)","v(HXC)"))
    else :
      logger.info("Found no hydrogen bonds to HD1")
    for typ in ["Donor: ","Acceptor: ","Any: ","Water: "] :
      if len(hid_hbonds[typ]) > 0 :
        ins = typ
        for item in hid_hbonds[typ] :
          logger.info("%10s\t%s%-5d%s %4s\t %8.2f %8.2f %8.2f %8.2f"%(ins,item[0].resname.capitalize(),item[0].residue,item[0].chain,item[0].name,item[1],item[2],item[3],item[4]))
          ins = ""
    
    logger.info("")
    if len(hie_hbonds["Donor: "]) > 0 or len(hie_hbonds["Acceptor: "]) > 0 or len(hie_hbonds["Any: "]) > 0 or len(hie_hbonds["Water: "]) > 0 :
      logger.info("**HE2**    \tResidue   Atom\t %8s %8s %8s %8s"%("r(X-H)","r(X-N)","v(NHX)","v(HXC)"))
    else :
      logger.info("Found no hydrogen bonds to HE2")
    for typ in ["Donor: ","Acceptor: ","Any: ","Water: "] :
      if len(hie_hbonds[typ]) > 0 :
        ins = typ
        for item in hie_hbonds[typ] :
          logger.info("%10s\t%s%-5d%s %4s\t %8.2f %8.2f %8.2f %8.2f"%(ins,item[0].resname.capitalize(),item[0].residue,item[0].chain,item[0].name,item[1],item[2],item[3],item[4]))
          ins = ""
    logger.info("-- -- -- --")

  raw_input("\nPress enter to continue")

  logger.info("")
  logger.info("Please select protonation state for the His residues (D=HID,E=HIE,P=HIP):")
  if user_sel == None :  
    hissel = []
    for i,h in enumerate(his) :
      sel = raw_input("His%d: "%h.serial)
      while not sel.upper() in ["D","E","P"] :
        print "Please enter D, E or P"
        sel = raw_input("His%d: "%h.serial)
      hissel.append(sel.upper())  
  else :
    hissel = user_sel

  logger.info("")
  logger.info("Selection: ")
  hlist = []
  for i,h in enumerate(his) :
    h.set_resname("HI%s"%hissel[i].upper())
    hlist.append("Hi%s%d"%(hissel[i].lower(),h.serial))
  logger.info(", ".join(hlist))

  return hissel


          
if __name__ == "__main__":

  # The only input argument is a PDB id
  pid = sys.argv[1].lower()

  # Setup the logger system
  logger = _setup_logger("%s_prep.log"%pid)

  # Try to read user selections from file or set to defaults
  sel = {}
  _read_sel(sel,pid)

  # If necessary download metadata and PDB-file from pdb.org
  pdbname = "%s.pdb"%pid
  metaname = "%s_meta.xml"%pid
  if sel["meta"] is None or (not os.path.isfile(metaname)) :
    _fetch("http://pdb.org/pdb/files/%s-noatom.xml"%pid.upper(),metaname,"meta data")
    sel["meta"] = metaname
  if sel["pdb"] is None or (not os.path.isfile(pdbname)):    
    _fetch("http://pdb.org/pdb/files/%s.pdb"%pid.upper(),pdbname,"PDB-file")
    sel["pdb"] = pdbname
  _write_sel(sel,pid)

  # Read PDB structure and obtain atoms,residues and chains
  pdbfile = PDBFile(filename="%s.pdb"%pid)
  
  # Perform various task that edits and describes the PDB file
  _new_section(wait=False)
  tasks = ["meta","buried","cyscmd","ignore_chains","ignore_atoms","ignore_het","his"]
  task_func = {"meta":_do_describe,"buried":_do_buried,"cyscmd":_do_cys,"ignore_chains":_do_chains,
               "ignore_het":_do_het,"ignore_atoms":_do_conformations,"his":_do_his}
  for task in tasks :
    sel[task] = task_func[task](pdbfile,pid,sel[task])
    _write_sel(sel,pid)
    _new_section()
  
  # Finally write out the prepared PDB file
  pdbfile.write("%s_prep.pdb"%pid,ter=True)
