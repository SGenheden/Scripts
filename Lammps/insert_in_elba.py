# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to insert an atomistic solute in box of ELBA molecules

Examples
--------
  insert_in_elba.py -o elba_toluene -b data.128dopc_4232wat -f forcefield.elba -s data.toluene -z {0..30} --respa
  insert_in_elba.py -b data.128dopc_4232wat_mod_whole -f forcefield.elba_mod -s data.kalp23_heavyh -o 128dopc_kalp23 --respa --resize -p kalp23_initial_placed.pdb
  insert_in_elba.py -b data.140popc_60chol_whole -f forcefield.140popc_60chol_mod2 -s data.b2_heavyh -o 140popc_60chol_b2 --respa --resize -p b2_initial_placed.pdb --dihfunc harmonic --dihfunc_box charmm --pairmix lj/charmm/coul/long/14:lj/charmm/coul/long
"""

import argparse
import sys
import os

import numpy as np

# Import the PDB and Lammps modules
thispath = os.path.dirname(os.path.abspath(__file__))
oneup = os.path.split(thispath)[0]
sys.path.insert(0,os.path.join(oneup,"Pdb"))
sys.path.insert(0,os.path.join(oneup,"Lammps"))
import lammps
import pdb

def _center_solute(solute_data,zcent) :

  # Center the solute in the membrane
  center = np.array([0.0,0.0,zcent])
  solutecenter = np.zeros(3)
  for atom in solute_data.atoms :
    solutecenter = solutecenter + atom.xyz
  diff = solutecenter / float(len(solute_data.atoms)) - center
  for atom in solute_data.atoms :
    atom.set_xyz(atom.xyz-diff)

if __name__ == '__main__' :

  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Program to insert an atomistic solute in an ELBA box")
  parser.add_argument('-o','--out',help="the name of the output ")
  parser.add_argument('-b','--box',help="the name of the box ")
  parser.add_argument('-f','--ff',help="name of the ELBA force field parameters",default="elba")
  parser.add_argument('-s','--solute',help="the name of the solute to be insert")
  parser.add_argument('-p','--pdb',help="a PDB-file of coordinates for the solute to insert")
  parser.add_argument('--respa',help="insert rRESPA hybrid flags into include file",action='store_true',default=False)
  parser.add_argument('--resize',help="resize box to fit solutes",action="store_true",default=False)
  parser.add_argument('-t','--type',choices=["cg2","cg1","aa"],help="the type of LAMMPS box, should be either 'cg2','cg1', 'aa'",default="cg1")
  parser.add_argument('-z','--zcent',nargs="+",type=float,help="the z-coordinate that the solute will be inserted at")
  parser.add_argument('--pairfunc',help="the pair function for the solute",default="lj/charmm/coul/long")
  parser.add_argument('--dihfunc',help="the dihedral function for the solute")
  parser.add_argument('--dihfunc_box',help="the dihedral function for the box")
  parser.add_argument('--pairmix',nargs="+",help="pair functions to use when mixing pairs")
  args = parser.parse_args()

  # Read solute data file
  solute_data = lammps.Datafile(filename=args.solute) 
  if args.pdb is not None :
    pdbfile = pdb.PDBFile(filename=args.pdb)
    for datom,patom in zip(solute_data.atoms,pdbfile.atoms) :
      datom.set_xyz(patom.xyz)

  # Read box data file and force field
  box_data = lammps.Datafile(filename=args.box)
  for atom in box_data.atoms :
    atom.ix = None
    atom.iy = None
    atom.iz = None
  box_ff = lammps.Includefile(filename=args.ff)

  if args.dihfunc_box is not None :
    for dih in box_ff.dihedralparams :
      if dih.func == "" : dih.func = args.dihfunc_box

  ELBA_FUNC = "lj/sf/dipole/sf"
  for pair in box_ff.pair_coeff :
    if pair.func == "" : pair.func = ELBA_FUNC
    if args.respa and pair.hybrid == -1 and pair.func == ELBA_FUNC: pair.hybrid = 2

  lj_hybrid = -1
  if args.respa : lj_hybrid = 1

  lj_func = {ELBA_FUNC:ELBA_FUNC}
  if args.respa :
    lj_hybrid = {ELBA_FUNC:1}
  else :
    lj_hybrid = None
  if args.pairmix is not None :
    for pm in args.pairmix :
      k,v = pm.split(":")
      lj_func[k] = v
      if args.respa : lj_hybrid[k] = -1

  # Extend the force field file with stuff from the datafile
  box_ff.extend_from_data(solute_data,lj_hybrid=-1,lj_func=args.pairfunc,lj_hybrid_mix=lj_hybrid,lj_func_mix=lj_func,ang_func="harmonic",dih_func=args.dihfunc)

  # Write a new force field file
  box_ff.write("forcefield."+args.out)

  # Combine the solute and box, put the solute at the end
  box_data.extend(solute_data)
  for atom in box_data.atoms :
    atom.kind = "cg/aa"

  if args.resize :
    xyz = np.zeros([len(box_data.atoms),3])
    for i,atom in enumerate(box_data.atoms) :
      xyz[i,:] = atom.xyz
    minxyz = xyz.min(axis=0)
    maxxzy  = xyz.max(axis=0)
  
    box_data.box[0] = minxyz[0]-2
    box_data.box[1] = minxyz[1]-2
    box_data.box[2] = minxyz[2]
    box_data.box[3] = maxxzy[0]-2
    box_data.box[4] = maxxzy[1]-2
    box_data.box[5] = maxxzy[2]

  if args.zcent is None :
    box_data.write("data."+args.out)
  else :
    for z in args.zcent :
      _center_solute(solute_data,z)
      box_data.write("data."+args.out+"_z%0.f"%z)      
  
