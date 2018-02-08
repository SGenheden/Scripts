# Author: Samuel Genheden samuel.genheden@gmail.com
"""
This script uses the prody library
"""

import argparse

from prody import *

def _calc_proj(trajectory, modes, filename, v1, v2):

    projection = calcProjection(trajectory, modes)
    with open(filename, "w") as f :
        f.write("#%.3f %.3f\n"%(v1,v2))
        for snap in projection :
            f.write("%.3f %.3f\n"%(snap[0], snap[1]))

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description="Program to project MD simulation on PCA modes")
    argparser.add_argument('-s','--struct',help="the filename of a PDB file")
    argparser.add_argument('-f','--file',help="the DCD trajectory")
    argparser.add_argument('-o','--out',help="output prefix")
    args = argparser.parse_args()

    structure = parsePDB(args.struct)
    trajectory = parseDCD(args.file)
    trajectory.setCoords(structure)
    trajectory.setAtoms(structure.calpha)
    trajectory.superpose()

    pca = PCA("PCA object")
    pca.buildCovariance(trajectory)
    pca.calcModes()

    v1 = calcFractVariance(pca[0])
    v2 = calcFractVariance(pca[1])

    n = int(0.5*trajectory.numConfs())
    trajsub1 = trajectory[:n]
    trajsub1.superpose()
    _calc_proj(trajsub1, pca[:2], args.out+"p1.out", v1, v2)

    trajsub2 = trajectory[n:]
    trajsub2.superpose()
    _calc_proj(trajsub2, pca[:2], args.out+"p2.out", v1, v2)
