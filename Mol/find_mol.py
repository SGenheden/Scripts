# Author: Samuel Genheden, samuel.genheden@gmail.com

import argparse
import os

from chemspipy import ChemSpider

from sgenlib import smiles

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Program to find molecules in ChemSpdie",)
    parser.add_argument('filename',help="the filename of a list")
    args = parser.parse_args()

    if os.getenv("SPIDERKEY") is None :
        print "SPIDERKEY environmental variable not set! Exit."
        quit()
    cs = ChemSpider(os.getenv("SPIDERKEY"))

    molecules = []
    with open(args.filename,"r") as f :
        molecules = [line.strip() for line in f.readlines()]

    for mol in molecules:
        hits = cs.search(mol)
        if len(hits) == 0 :
            print mol+"\t!!"
        else :
            """try :
                print "//".join([h.common_name for h in hits])
            except :
                print mol+"\t!!!"
                """
            print "%d"%hits[0].csid+"\t"+hits[0].common_name
