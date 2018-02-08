# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to calculate RDF in a two-component system
"""

import argparse

import MDAnalysis
from MDAnalysis.analysis import rdf

def _anal_rdf(u, sel1, sel2, exl, filename):

    g1 = u.select_atoms(sel1)
    g2 = u.select_atoms(sel2)
    rdfo = rdf.InterRDF(g1, g2, exclusion_block=exl, step=50)
    rdfo.run()
    with open(filename, "w") as f:
        for b, r in zip(rdfo.bins, rdfo.rdf):
            f.write("%.3f %.3f\n"%(b,r))

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description="Script to count hydrogen bonds")
    argparser.add_argument('-f', '--file', help="the trajectory file.")
    argparser.add_argument('-s', '--struct', help="a structure file")
    argparser.add_argument('--nwat', type=int, help="the number of water atoms", default=4)
    argparser.add_argument('-out','--out', help="the output file",  default="rdf")

    args = argparser.parse_args()

    u = MDAnalysis.Universe(args.struct, args.file)
    #_anal_rdf(u, "resname 0GA or resname 1GA", "resname SOL", None, args.out+"_tw.txt")
    #_anal_rdf(u, "resname 0GA or resname 1GA", "resname 0GA or resname 1GA", (45, 45),  args.out+"_tt.txt")
    _anal_rdf(u, "resname SOL", "resname SOL", (args.nwat, args.nwat), args.out+"_ww.txt")

    #_anal_rdf(u, "resname SOL and name O*", "resname SOL and name O*", (1,1), args.out+"_ww-o.txt")
    #_anal_rdf(u, "resname SOL and name H*", "resname SOL and name H*", (2,2), args.out+"_ww-h.txt")
    #_anal_rdf(u, "resname SOL and name O*", "resname SOL and name H*", None, args.out+"_ww-oh.txt")

    #_anal_rdf(u, "(resname 0GA or resname 1GA) and name O*", "(resname 0GA or resname 1GA) and name O*", (11, 11),  args.out+"_tt-o.txt")
    #_anal_rdf(u, "(resname 0GA or resname 1GA) and name H*", "(resname 0GA or resname 1GA) and name O*", (22, 11),  args.out+"_tt-oh.txt")
    #_anal_rdf(u, "(resname 0GA or resname 1GA) and name H*", "(resname 0GA or resname 1GA) and name H*", (22, 22),  args.out+"_tt-hh.txt")
    #_anal_rdf(u, "name H2O or name H3O or name H4O or name H6O", "name H2O or name H3O or name H4O or name H6O", (8, 8), args.out+"_tt-hhx.txt")

    #_anal_rdf(u, "(resname 0GA or resname 1GA) and name O*", "resname SOL and name H*", None, args.out+"_tw-oh.txt")
    #_anal_rdf(u, "(resname 0GA or resname 1GA) and name H*", "resname SOL and name O*", None, args.out+"_tw-ho.txt")

    #_anal_rdf(u, "name O2 or name O3 or name O4", "resname SOL and name H*", None, args.out+"_tw-ohw.txt")
    #_anal_rdf(u, "name O6", "resname SOL and name H*", None, args.out+"_tw-o6hw.txt")
    #_anal_rdf(u, "name O1", "resname SOL and name H*", None, args.out+"_tw-o1hw.txt")
    #_anal_rdf(u, "name O5", "resname SOL and name H*", None, args.out+"_tw-o5hw.txt")
    #_anal_rdf(u, "name O2 or name O3 or name O4", "resname SOL and name O*", None, args.out+"_tw-oow.txt")
    #_anal_rdf(u, "name O6", "resname SOL and name O*", None, args.out+"_tw-o6ow.txt")
    #_anal_rdf(u, "name O1", "resname SOL and name O*", None, args.out+"_tw-o1Ow.txt")
    #_anal_rdf(u, "name O5", "resname SOL and name O*", None, args.out+"_tw-o5Ow.txt")

    #_anal_rdf(u, "name H2O or name H3O or name H4O or name H6O", "resname SOL and name O*", None, args.out+"_tw-oxow.txt")
