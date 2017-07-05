import argparse
import os
import sys

from sgenlib import pdb
from sgenlib import gmx

if __name__ == "__main__":

  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Program sort a structure file according to topology")
  parser.add_argument('-f','--file',help="the input gro-file")
  parser.add_argument('-p','--topol',help="the topology file")
  parser.add_argument('-o','--out',help="the output gro-file",default="reordered.gro")
  args = parser.parse_args()

  ref = gmx.TopFile(args.topol)
  pdbfile = pdb.PDBFile(args.file)

  resnames = []
  atomnames = []
  reslen = {}
  for mol in ref.mollist :
    for moltype in ref.moleculetypes :
      if mol.lower() != moltype.name.lower() : continue
      found = True
      refatoms = moltype.atoms
      break
    if not found : continue
    n = ref.molecules[mol]
    resnames.extend([refatoms[0].resname]*n)
    atomnames.extend([[a.name for a in refatoms] for i in range(n)])
    reslen[refatoms[0].resname.strip()] = len(refatoms)

  print reslen
  for i,res in enumerate(pdbfile.residues) :
    if res.resname.strip() not in reslen :
      raise Exception("Found %s in structure but not in topology"%res.resname.strip())
    else :
      slen = len(res.atoms)
      tlen = reslen[res.resname.strip()]
      if slen != tlen :
        if slen % tlen != 0 :
          raise Exception("Found a residue with incompatible length (%s, %s, %s)"%(tlen,slen,res.resname))
        else :
          pdbfile.split_residue(i,slen/tlen)

  pdbfile.reorder(resnames,atomnames)
  pdbfile.write_gro(args.out)

  """
  print "Reordered"
  print "%5d"%len(pdbfile.atoms)

  natom = 0
  for mol in ref.mollist :
    for ri,residue in enumerate(pdbfile.residues) :
      found = False
      for moltype in ref.moleculetypes :
        if mol.lower() != moltype.name.lower() : continue
        if moltype.name.find(residue.resname) != -1 or moltype.atoms[0].resname == residue.resname :
          found = True
          refatoms = moltype.atoms
          break
      if not found : continue
      for ratom in refatoms :
        found = False
        for atom in residue.atoms :
          if ratom.name.strip() == atom.name.strip() :
            found = True
            natom = natom + 1
            if natom > 99999 : natom = natom - 99999
            print "%5d%-5s%5s%5d%8.3f%8.3f%8.3f"%((ri+1 if ri+1<=9999 else ri+1-9999),atom.resname,atom.name,natom,atom.x/10.0,atom.y/10.0,atom.z/10.0)
        if not found :
          print "Fatal, could not find topology atom %s in structure file"%(ratom.name)
          quit()
  print "%8.3f%8.3f%8.3f"%(pdbfile.box[0]/10.0,pdbfile.box[1]/10.0,pdbfile.box[2]/10.0)
  """
