# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to make input to the colvars library to perform umbrella sampling simulations
with Lammps

It will produce an output file for each value of the --zdepth argument

Example:
  make_colvars.py data.elba_toluene_z0 -m 1 128 -s 129 -z {0..30}
"""

import argparse

from sgenlib import lammps

if __name__ == '__main__' :

  # Command-line input
  parser = argparse.ArgumentParser(description="Making colvars input for an umbrella simulation")
  parser.add_argument('file',help="a template data file.",default="")
  parser.add_argument('-o','--out',help="the output file prefix.",default="colvars.")
  parser.add_argument('-m','--membrane',type=int,nargs=2,help="the first and last molecule of the membrane")
  parser.add_argument('-s','--solute',type=int,help="the molecule id of the solute")
  parser.add_argument('-z','--zdepth',nargs="+",type=float,help="the z-depth")
  parser.add_argument('-w','--weight',type=float,help="the weight of the umbrella restraint",default=2.5)
  args = parser.parse_args()

  strtempl = """
colvarsTrajFrequency 150
colvarsRestartFrequency 0
 
colvar {
  name solute
   
  width 1.0

  distanceZ {
    main {
      # Solute atoms
      atomNumbersRange %d-%d
    }
    ref { 
      # Lipid atoms     
      atomNumbersRange %d-%d
      #dummyAtom (0.0,0.0,0.0)
    }
    axis (0,0,1)
  }
}

harmonic {
  colvars solute
  forceConstant %.3f
  centers %.1f
  outputCenters on
} 

"""

  
  # Find the first and last atom of the membrane and solute
  mfirst = 10000
  mlast  = -10000
  sfirst = 10000
  slast  = -10000
  for atom in lammps.Datafile(args.file).atoms :
    if atom.molecule >= args.membrane[0] and atom.molecule <= args.membrane[1] :
      mfirst = min(atom.idx,mfirst)
      mlast = max(atom.idx,mlast)
    elif atom.molecule == args.solute :
      sfirst = min(atom.idx,sfirst)
      slast = max(atom.idx,slast)

  # Create an output file for each depth
  for z in args.zdepth :
    with open("%sz%0.f"%(args.out,z),"w") as f:
      f.write(strtempl%(sfirst,slast,mfirst,mlast,args.weight,z))
          
