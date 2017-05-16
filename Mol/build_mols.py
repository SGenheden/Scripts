# Author: Samuel Genheden, samuel.genheden@gmail.com

"""
Program to build molecules by looking up there SMILES string in
ChemSpider and then converting the SMILES to an xyz file.

The SPIDERKEY environmental variable needs to be set ChemSpider access key.

If the molecule could not be found in ChemSpider the molecule name is
printed with a double exclamation mark (!!) appended.

If the SMILES could not be converted successfully the molecule name is
printed with a single exclamation mark (!) appended.

If multiple hits are found in ChemSpider the first one is converted
to a structure

Examples
--------
  build_mols.py molecules
"""

import argparse
import os

from chemspipy import ChemSpider

from sgenlib import smiles

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Program to build molecules",)
    parser.add_argument('filename',help="the filename of a list")
    args = parser.parse_args()

    if os.getenv("SPIDERKEY") is None :
        print "SPIDERKEY environmental variable not set! Exit."
        quit()
    cs = ChemSpider(os.getenv("SPIDERKEY"))

    lines = []
    with open(args.filename,"r") as f :
        lines = [line.strip() for line in f.readlines()]
    molecules = sorted(list(set(lines)),cmp=lambda x,y: cmp(lines.index(x),lines.index(y)))

    for mol in molecules:
        hits = cs.search(mol)
        if len(hits) == 0 :
            print mol+"\t!!"
        else :
            molsmiles = hits[0].smiles
            mol2 = mol.strip()
            mol2.replace(" ","_")
            if smiles.convert2xyz(molsmiles,mol2+".xyz",verbose=False) is not None :
                print mol+"\t"+molsmiles
            else:
                print mol+"\t!"
