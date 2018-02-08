# Author: Samuel Genheden samuel.genheden@gmail.com
"""
Calculate interaction energy between the protein and a ligand, by extracing
Gromacs xvg-files
"""

import argparse

import numpy as np

from sgenlib import pdb
from sgenlib import parsing

def _make_stats(data) :

    mean = data.mean()
    std = data.std()/np.sqrt(data.shape[0])
    nhalf = int(0.5*data.shape[0])
    drift = np.abs(data[:nhalf].mean() - data[nhalf:].mean())
    return "%.3f\t%.3f\t%.3f"%(mean, std, drift),[mean, std, drift]

if __name__ == '__main__' :

  # Command-line input
  parser = argparse.ArgumentParser(description="Calculate ligand-protein interaction energies")
  parser.add_argument('files',nargs="+",help="the Gromacs xvg-files")
  parser.add_argument('-o','--out',help="the output file")
  parser.add_argument('-p','--pdb',help="the protein",default=[])
  parser.add_argument('-r','--repeats', help="repeat replacement", default=["r1_", "r2_", "r3_", "r4_", "r5_", "r6_", "r8_", "r9_", "r10_"])
  args = parser.parse_args()

residues = ["CU", "CL"]
for res in pdb.PDBFile(args.pdb).residues[:235] :
    residues.append("%s%d"%(res.resname.capitalize(),res.serial))

all_data = None
for ri, replacement in enumerate(args.repeats):
    r_data = []
    if ri == 0 :
        outname = args.out
    else:
        outname = args.out.replace(args.repeats[0], replacement)
    with open(outname, "w") as f :
        if ri == 0 :
            files = args.files
        else :
            files = [file.replace(args.repeats[0], replacement) for file in args.files]
        for file in files :
            data = parsing.parse2ndarray(file)
            (nrows, nitems) = data.shape
            for i in range(1, nitems, 2):
                ele_str, ele = _make_stats(data[:, i])
                vdw_str, vdw = _make_stats(data[:, i+1])
                tot_str, tot = _make_stats(data[:, i] + data[:, i+1])

                ele.extend(vdw)
                ele.extend(tot)
                r_data.append(ele)

                f.write("%s\t%s\t%s\t%s\n"%(residues[i-1], ele_str, vdw_str, tot_str))

        if ri == 0 :
            r_data = np.asarray(r_data)
            nrow, ncol = r_data.shape
            all_data = np.zeros([len(args.repeats), nrow, ncol])
            all_data[0, :, : ] = r_data
        else :
            all_data[ri, :, : ] = np.asarray(r_data)

with open(args.out.replace(args.repeats[0], "all_"), "w") as f :
    mean_data = all_data.mean(axis=0)
    for i in range(0, 7, 3):
        mean_data[:, i+1] = np.std(all_data[:, i])/np.sqrt(len(args.repeats))
    for i, residue in enumerate(residues) :
        f.write("%s\t%s\n"%(residue,"\t".join("%.3f"%d for d in mean_data[i, :])))
