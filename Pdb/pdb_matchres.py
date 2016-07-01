# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to match the residue numbers of one structure
with another.
"""

import sys
import argparse

from sgenlib import pdb

if __name__ == "__main__":


  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Program to match residue number of 2 PDB files")
  parser.add_argument('-reference','--reference',help="the reference PDB file")
  parser.add_argument('-mobile','--mobile',help="the mobile PDB file")
  parser.add_argument('-ranges','--ranges',nargs="+",help="list of residues ranges to match")
  parser.add_argument('-offset','--offset',default=0,type=int,help="reference offset")
  args = parser.parse_args()

  pdb1 = pdb.PDBFile(args.reference, renumber=False)
  pdb2 = pdb.PDBFile(args.mobile, renumber=False)

  naa1 = sum([1  for r in pdb1.residues if r.resname in pdb.std_aa_names])
  naa2 = sum([1  for r in pdb2.residues if r.resname in pdb.std_aa_names])
  if naa1 != naa2 :
    raise Exception("The two PDB files does not have the same number of residues (%d, %d)."%(naa1,naa2))

  if args.offset != 0 :
    for res in pdb1.residues : res.serial += args.offset

  print "[Residue match]"
  print "Command= "+" ".join(sys.argv)
  print "Numbers-ref = %s"%" ".join(["%d"%(res.serial) for res in pdb1.residues])
  print "Numbers-mob = %s"%" ".join(["%d"%res.serial for res in pdb2.residues if res.resname in pdb.std_aa_names])
  match = {}
  for r1,r2 in zip(pdb1.residues,pdb2.residues) :
    print "%d = %d"%(r1.serial,r2.serial)
    match[r1.serial] = r2

  print "[Residue names]"
  print "Names = %s"%" ".join(res.resname for res in pdb1.residues)
  print "Codes = %s"%" ".join(pdb.codes[res.resname.lower()] for res in pdb1.residues)

  if args.ranges is not None :
    print "[Residue ranges]"
    print "Range = %s"%" ".join(args.ranges)
    deflist = []
    defalist = []
    for r in args.ranges :
      r1,r2 = map(int,r.split("-"))
      res = (match[r1].serial,match[r2].serial)
      atm = (match[r1].atoms[0].serial,match[r2].atoms[-1].serial)
      print "%s = %s %s"%(r,"%s-%s"%res,"%s-%s"%atm)
      deflist.append("[%s, %s]"%res)
      defalist.append("[%s, %s]"%atm)
    print "Def-list = [%s]"%(",".join(deflist))
    print "Def-atom-list = [%s]"%(",".join(defalist))
