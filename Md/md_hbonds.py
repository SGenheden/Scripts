# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to count hydrogen bonds in a protein

The analysis is done on all atoms but also on only backbone
and only side chain atoms

Examples
--------
md_hbonds.py -f md2.xtc -s md1.gro
"""

import argparse

import MDAnalysis
from MDAnalysis.analysis import hbonds

class HydrogenBondAnalysis_mainchain(hbonds.HydrogenBondAnalysis):
    DEFAULT_DONORS = {
        'CHARMM27': tuple(set([
            'N']))}
    DEFAULT_ACCEPTORS = {
        'CHARMM27': tuple(set([
            'O']))}

class HydrogenBondAnalysis_sidechain(hbonds.HydrogenBondAnalysis):
    DEFAULT_DONORS = {
        'CHARMM27': tuple(set([
            'NE', 'NH1', 'NH2', 'ND2', 'SG', 'NE2', 'ND1', 'NZ', 'OG', 'OG1', 'NE1', 'OH']))}
    DEFAULT_ACCEPTORS = {
        'CHARMM27': tuple(set([
            'OD1', 'OD2', 'SG', 'OE1', 'OE1', 'OE2', 'ND1', 'NE2', 'SD', 'OG', 'OG1', 'OH']))}

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description="Script to count hydrogen bonds")
    argparser.add_argument('-f', '--file', help="the trajectory file.")
    argparser.add_argument('-s', '--struct', help="a structure file")
    args = argparser.parse_args()

    u=MDAnalysis.Universe(args.struct,args.file)
    h_all = hbonds.HydrogenBondAnalysis(u,"protein","protein",
                update_selection1=False, update_selection2=False)
    results_all = h_all.run()
    cnt = h_all.count_by_time().count
    print "All:\t%.3f\t%.3f"%(cnt.mean(),cnt.std())

    h_main = HydrogenBondAnalysis_mainchain(u,"protein","protein",
                update_selection1=False, update_selection2=False)
    results_main = h_main.run()
    cnt = h_main.count_by_time().count
    print "Main:\t%.3f\t%.3f"%(cnt.mean(),cnt.std())

    h_side = HydrogenBondAnalysis_sidechain(u,"protein","protein",
                update_selection1=False, update_selection2=False)
    results_side = h_side.run()
    cnt = h_side.count_by_time().count
    print "Side:\t%.3f\t%.3f"%(cnt.mean(),cnt.std())
