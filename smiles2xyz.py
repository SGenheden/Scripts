# Author: Samuel Genheden, samuel.genheden@gmail.com

"""
Program to convert a SMILES string to 3D coordinates using
one of two web services. 

Examples
--------
  smiles2xyz.py COH
  smiles2xyz.py -o ethanol.cyx CCOH
"""

import sys
import urllib
import os
import argparse
 
def cactus(smiles) :
  """
  Convert SMILES using Cactus service
  
  Parameters
  ----------
  smiles : string
    the SMILES string
  
  Returns
  -------
  string
    the path to the downloaded SDF file or none if conversion failed
  """
  site = "cactus.nci.nih.gov"
  try:
    page = urllib.urlopen("http://%s/cgi-bin/translate.tcl?smiles=%s&format=sdf&astyle=kekule&dim=3D&file="%(site, smiles))
  except : 
    return None
  for line in page:
    if "Click here" in line  and 'a href="' in line :
	  dummy1, url, dummy2 = line.split('"')
	  try :
	    path, header = urllib.urlretrieve("http://%s%s"%(site,url))
	  except :
	    return None
	  return path

def indiana(smiles) :
  """
  Convert SMILES using Indiana service
  
  Parameters
  ----------
  smiles : string
    the SMILES string
  
  Returns
  -------
  string
    the path to the downloaded SDF file or none if conversion failed
  """
  try :
    path, header = urllib.urlretrieve("http://cheminfov.informatics.indiana.edu/rest/thread/d3.py/SMILES/%s" % smiles)
  except :
    return None
  return path
  
def sdf2xyz(sdffile,xyzfile) :
  """
  Read a SDF file and convert it to xyz-format
  
  Parameters
  ----------
  sdffile : string
    the name of the SDF file
  xyzfile : string
    the name of the XYZ-file
  """
  lines = open(sdffile,"r").readlines()
  if len(lines) == 0 : return
  with open(xyzfile,"w") as f :
    natom = int(lines[3].strip().split()[0])
    f.write("%d\n0\n"%natom)
    for line in lines[4:4+natom] :
      cols = line.strip().split()
      f.write("%s %s %s %s\n"%(cols[3],cols[0],cols[1],cols[2]))
   
if __name__ == "__main__":

  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Program to convert a SMILES string to 3D-coordinates in xyz format",)
  parser.add_argument('smiles',help="the SMILES string")
  parser.add_argument('-o','--out',help="the output xyz-file, default='mol.xyz'",default="mol.xyz")
  parser.add_argument('-w','--web',choices=["cactus","indiana","both"],help="the name of the web service used to convert the SMILES, should be either 'cactus','indiana', or 'both'",default="both")
  args = parser.parse_args()

  # Try to convert the SMILES string using one of two web services
  service = ""
  if args.web == "cactus" :
    path = cactus(args.smiles.strip())
    if len(open(path,"r").readlines()) == 0 : path = False
    service = "cactus"
  elif args.web == "indiana" :
    path = indiana(args.smiles.strip())
    if len(open(path,"r").readlines()) == 0 : path = False
    serivce = "indiana"
  else :
    path = indiana(args.smiles.strip())
    if len(open(path,"r").readlines()) == 0 : path = False
    service = "indiana"
    if not path :
      path = cactus(args.smiles.strip())
      service = "cactus"
  
  # If the conversion was successful convert it to xyz-format 
  text = {"indiana":"the smi23d web service provided by the Chemical Informatics and Cyberinfrastructure Collaboratory at Indiana University","cactus":" the SMILES translator provided by the National Cancer Institute CADD group"}  
  link = {"indiana":"http://www.chembiogrid.info/projects/proj_ws_all.html","cactus":"http://cactus.nci.nih.gov/translate/"}
  if path :
    sdf2xyz(path,args.out)
    print "SMILES converted sucessfully using %s"%text[service]
    print "Web page: %s"%link[service]
  else :
    print "Unable to convert SMILES."
