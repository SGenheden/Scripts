# Author: Samuel Genheden, samuel.genheden@gmail.com

"""
Program to draw 2D structures on a grid with RDKit
"""

import argparse

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

if __name__ == '__main__' :

    parser = argparse.ArgumentParser(description="Draw molecules")
    parser.add_argument('files',nargs="+",help="the data files")
    parser.add_argument('-o','--out',help="the output prefix",default="mols")
    args = parser.parse_args()

    mols = [Chem.SDMolSupplier(filename)[0] for filename in args.files]
    for m in mols: tmp=AllChem.Compute2DCoords(m)
    for i in range(0,len(mols),30):
        img=Draw.MolsToGridImage(mols[i:i+30],molsPerRow=5,subImgSize=(200,200),legends=args.files[i:i+30])
        img.save("%s%d.png"%(args.out,i+1))
