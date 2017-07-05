# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to scale CG-AA pair interactions in LAMMPS inclusion file

The CG-AA pairs need to have a comment including "AA-CG mixed" to be regonized
as such.

The --beads argument is a list of bead specification, each one contains a comma-separated
list of atom types. The pairs between these atom types and AA types will be scaled.

The --scale argument is the scaling factor and --param can be either e, s, or c,
indicating scaling epsilon, sigma (vdW parameters) or Coulomb potential.

The default output is the input include file with a "_scaled" string appended

Example:
 (standard scaling in hybrid membrane simulation)
  scale_elba_pairs.py forcefield.elba_toluene -b 2,3 4,5 -s 0.5 1.75 -p e c
  (standard scaling in hybrid liquid simulation)
  scale_elba_pairs.py forcefield.elba_toluene -b 6 -s 0.9 -p e
"""

import argparse

import numpy as np

from sgenlib import lammps

if __name__ == '__main__' :

  # Command-line input
  parser = argparse.ArgumentParser(description="Scaling CG-AA interaction in LAMMPS inclusion file")
  parser.add_argument('file',help="the inclusion file.")
  parser.add_argument('-o','--out',help="the output file.")
  parser.add_argument('-b','--beads',nargs="+",help="the CG atom type groups",default=[])
  parser.add_argument('-s','--scale',type=float,nargs="+",help="the scaling factors",default=[])
  parser.add_argument('-p','--param',choices=["e","s","c"],nargs="+",help="the parameter to scale, either e(psilon), s(igma), or c(oulomb)",default=[])
  args = parser.parse_args()

  if args.file is None :
    print "No input file specified. Exiting!"
    quit()

  infile = lammps.Includefile(filename=args.file)

  # Find all mixed pairs
  mixed_pairs = []
  for pair in infile.pair_coeff :
    if (pair.comment.find("AA-CG mixed") > -1 or pair.comment.find("Mixed using Lorentz-Berthelot rules") > -1 )  :
      mixed_pairs.append(pair)
  lj_same = infile.lj_same_dict()

  # Loop over all the bead specification and change the given parameter
  for beadspec,fac,param in zip(args.beads,args.scale,args.param) :
    beads = [int(b) for b in beadspec.strip().split(",")]
    for pair in mixed_pairs :
      if pair.iatom in beads :
        if param.lower()[0] == "e" :
          new = fac*np.sqrt(lj_same[pair.iatom].epsilon*lj_same[pair.jatom].epsilon)
          pair.epsilon = new
        elif param.lower()[0] == "s" :
          new = fac*0.5*(lj_same[pair.iatom].sigma+lj_same[pair.jatom].sigma)
          pair.sigma = new
        elif param.lower()[0] == "c" :
          pair.scale = fac

  # Write output
  if args.out is None :
    args.out = args.file+"_scaled"
  infile.write(args.out)
