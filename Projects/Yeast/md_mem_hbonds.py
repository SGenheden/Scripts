# Author: Samuel Genheden samuel.genheden@gmail.com
"""
Script to count hydrogen bonds
"""

import argparse

import MDAnalysis
from MDAnalysis.analysis import hbonds

class HydrogenBondAnalysis_lipids(hbonds.HydrogenBondAnalysis):
    DEFAULT_DONORS = {
        'CHARMM27': tuple(set([
            'O2F', 'O4S', 'O1S', 'O2', 'O3', 'O4', 'O5', 'O6', 'NF', 'OH2']))}
    DEFAULT_ACCEPTORS = {
        'CHARMM27': tuple(set([
            'OF', 'O2F', 'O4S', 'O1S', 'O13', 'O14', 'O2', 'O3', 'O4', 'O5', 'O6', 'O22', 'O32', 'OH2', 'O1']))}

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description="Script to count hydrogen bonds")
    argparser.add_argument('-f', '--file', help="the trajectory file.")
    argparser.add_argument('-s', '--struct', help="a structure file")
    argparser.add_argument('--sel1', help="the first selection")
    argparser.add_argument('--sel2', help="the second selection")
    argparser.add_argument('-out','--out', help="the output pickle",  default="mem_hbonds.pickle")

    args = argparser.parse_args()

    u=MDAnalysis.Universe(args.struct,args.file)
    h = HydrogenBondAnalysis_lipids(u,args.sel1,args.sel2,
                update_selection1=False, update_selection2=False)
    results = h.run()
    h.generate_table()
    h.save_table(args.out)
