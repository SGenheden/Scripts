# Author: Samuel Genheden, samuel.genheden@gmail.com

"""
Program to convert a GRO file to PDB and in the process add CONECT records
"""

import sys,os
from sgenlib import pdb
import argparse

def make_con(filename) :

  lines = open(filename,"r").readlines()
  con_raw = []
  all_ind = []
  for line in lines :
    cols = line.strip().split()
    con_raw.append((int(cols[0]),int(cols[1])))
    all_ind.append(int(cols[0]))
    all_ind.append(int(cols[1]))
  natom = max(all_ind)
  con = [[] for i in range(natom)]
  for c in con_raw :
    con[c[0]-1].append(c[1]-1)
    con[c[1]-1].append(c[0]-1)
  return con

def add_con(pdbfile,f) :

  if prot_con :
    for i,c in enumerate(prot_con) :
      c2 = set(c)
      f.write("CONECT%5d%s\n"%(i+1,"".join(["%5d"%(j+1) for j in c2])))
  if het_con :
   for residue in pdbfile.residues :
     for resname in het_con.keys() :
       if residue.resname[:3] == resname[:3] :
         for i,c in enumerate(het_con[resname]) :
           c2 = set(c)
           f.write("CONECT%5d%s\n"%(residue.atoms[0].serial+i,"".join(["%5d"%(residue.atoms[0].serial+j) for j in c2])))


if __name__ == '__main__' :

  # Command-line input

  parser = argparse.ArgumentParser(description="Converting GRO files to PDB and add connectivity")
  parser.add_argument('-f','--files',nargs="+",help="the input GRO-files.",default=[])
  parser.add_argument('-pb','--pbonds',help="the connectivity information for the protein",default="")
  parser.add_argument('-hb','--hbonds',nargs="+",help="the connectivity information for the het residues",default=[])
  args = parser.parse_args()

  if len(args.files) == 0 :
    print "Nothing to do. Exiting."
    quit()

  prot_con = None
  if args.pbonds == "" :
    print "No connectivity information will be added for the protein"
  else :
    prot_con = make_con(args.pbonds)

  het_con = None
  if len(args.hbonds) == 0 :
    print "No connectivity information will be added for the hetero residues"
  else :
    het_con = {}
    for filename in args.hbonds :
      filebase,fileext = os.path.splitext(filename)
      het_con[filebase.upper()] = make_con(filename)

  print het_con

  for filename in args.files :
    filebase,fileext = os.path.splitext(filename)
    pdbfile = pdb.PDBFile()
    pdbfile.read(filename,gro=True)
    print len(pdbfile.chains)
    pdbfile.write(filebase+".pdb",ter=True,add_extra=add_con)
