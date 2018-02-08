# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to count hydrogen bonds in a two-component system

Examples:
    python calc_hbonds_twocomp.py -f r1_md2_whole.xtc -s r1_md2.tpr --sel "all" -out r1_hbonds.pickle
    python calc_hbonds_twocomp.py -f template.xyz wt1kalpotAllAtoms.xyz -s template.pdb --angle 150
"""

import argparse

import MDAnalysis
from MDAnalysis.analysis import hbonds

class HydrogenBondAnalysis_twocomp(hbonds.HydrogenBondAnalysis):
    DEFAULT_DONORS = {
        'CHARMM27': tuple(set([
            'OW','O2','O3','O4','O6']))}
    DEFAULT_ACCEPTORS = {
        'CHARMM27': tuple(set([
            'OW','O2','O3','O4','O6']))}


if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description="Script to count hydrogen bonds")
    argparser.add_argument('-f', '--file', nargs="+", help="the trajectory file.")
    argparser.add_argument('-s', '--struct', help="a structure file")
    argparser.add_argument('--sel', help="the first selection", default="all")
    argparser.add_argument('-out','--out', help="the output pickle",  default="hbonds.pickle")
    argparser.add_argument('--distance', type=float, help="the cut-off distance",  default=3)
    argparser.add_argument('--angle', type=float, help="the cut-off angle",  default=120)

    args = argparser.parse_args()
    #MDAnalysis.start_logging()

    u=MDAnalysis.Universe(args.struct,*args.file)
    h = HydrogenBondAnalysis_twocomp(u,args.sel,args.sel, distance=args.distance, angle=args.angle,
                update_selection1=False, update_selection2=False) #, start=200, step=25)
    results = h.run()
    h.generate_table()
    h.save_table(args.out)
