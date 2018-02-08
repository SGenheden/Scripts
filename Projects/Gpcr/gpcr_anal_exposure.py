# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to draw roughness (fractal) pies on structures

No arguments are necessary, all structures are taken from standard locations
"""

import argparse
import os
import sys

import numpy as np
import scipy.stats as stats

import matplotlib
if not "DISPLAY" in os.environ or os.environ["DISPLAY"] == "" :
  matplotlib.use('Agg')
import matplotlib.pylab as plt

import  gpcr_lib
# Import the calc_surf program
thispath = os.path.dirname(os.path.abspath(__file__))
oneup = os.path.split(thispath)[0]
sys.path.insert(0,oneup)
import calc_surf

if __name__ == '__main__' :

    # Command-line input
    parser = argparse.ArgumentParser(description="Analysing residue exposure")
    parser.add_argument('-f','--folder', help="the  folder with the residue contacts")
    parser.add_argument('-p','--probe',type=float,help="the probe size",default=2.4)
    parser.add_argument('-c','--percentile',type=float,help="the percentile to consider",default=50)
    args = parser.parse_args()

    mols = "b2 b2_a a2a a2a_a".split()

    data = {}
    for mol in mols:
        xray, aastruct = gpcr_lib.load_xray(mol, loadsigma=True, loadaa=True)
        radii = np.asarray([calc_surf.bornradii[atom.element().upper()] for atom in aastruct.atoms])
        xyzrname = calc_surf.write_xyzr(aastruct.xyz,radii)
        surf = calc_surf.calc_surf(xyzrname, aastruct.xyz.shape[0], args.probe)
        res, contactprob = gpcr_lib.read_rescontacts(args.folder, mol, percentile=args.percentile, returnprobs=True)
        print len(res)
        mdata = {}
        for tres, s in zip(xray.template.residues, surf):
            if tres in res :
                mdata[tres] = (s[0], contactprob[res==tres][0])
        data[mol] = mdata
        with open("ses_%s"%mol,"w") as f :
            for s, r in zip(surf,xray.template.residues):
                f.write("%d %.3f\n"%(r, s[0]))

    diff = []
    for res in data["b2"]:
        if res in data["b2_a"]:
            diff.append([data["b2_a"][res][0]-data["b2"][res][0],
                    data["b2_a"][res][1]-data["b2"][res][1]])
    diff = np.asarray(diff)
    print "B2: %.3f %.3f"%(np.corrcoef(diff[:,0],diff[:,1])[1,0],
        stats.kendalltau(diff[:,0],diff[:,1])[0])


    diff = []
    for res in data["a2a"]:
        if res in data["a2a_a"]:
            diff.append([data["a2a_a"][res][0]-data["a2a"][res][0],
                    data["a2a_a"][res][1]-data["a2a"][res][1]])
    diff = np.asarray(diff)
    print "A2a: %.3f %.3f"%(np.corrcoef(diff[:,0],diff[:,1])[1,0],
        stats.kendalltau(diff[:,0],diff[:,1])[0])
