# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to solvate a MARTINI CG or GAFF/CG system

Assumes last include is the solute itp file

Examples:
    solvate_martini.py -b octanol_box.gro -s eth_cg_box.gro -os eth_solvated.gro -i martini_v2.2.itp martini_v2.0_solvents.itp eth_cg.itp -m OCO -ot eth_system.top
    solvate_martini.py -b water_box.gro -s eth_cg_box.gro -os eth_solvated.gro -i martini_v2.2.itp eth_cg.itp -m W -ot eth_system.top

"""

import argparse
import os

from sgenlib import pdb
from sgenlib import ambertools

topfile_str="""

[ System ]
Solvated %s

[ Molecules ]
%s     1
%s     %d
"""

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Program solvate MARTINI solute")
    parser.add_argument('-i','--includes',nargs="+", help="the include files", default=[])
    parser.add_argument('-b','--box', help="the solvent box file")
    parser.add_argument('-s','--solute', help="the solute file")
    parser.add_argument('-m','--molecule', help="the molecule name of the solvent")
    parser.add_argument('-ot','--outtop',help="the output top-file",default="system.top")
    parser.add_argument('-os','--outstruct',help="the output gro-file",default="solvated.gro")
    parser.add_argument('--makendx',action="store_true",help="write a simple ndx file",default=False)
    args = parser.parse_args()

    solvate_cmd = "$GMX/gmx solvate -cp %s -cs %s -o %s >& solv.log"%(
                    args.solute, args.box, args.outstruct)
    ambertools.run_program("cmd", solvate_cmd)

    solute = pdb.PDBFile(args.solute)
    solvbox = pdb.PDBFile(args.box)
    solvated = pdb.PDBFile(args.outstruct)
    natomperres = len(solvbox.atoms) / len(solvbox.residues)
    nsolvent = (len(solvated.atoms) - len(solute.atoms)) / natomperres

    with open(args.includes[-1], "r") as f :
        line = f.readline()
        while line and line.find("moleculetype") == -1 :
            line = f.readline()
        line = f.readline()
        line = f.readline()
        solutename, dummy = line.strip().split()

    with open(args.outtop, "w") as f :
        for i in args.includes :
            f.write('#include "%s"\n'%i)
        f.write(topfile_str%(os.path.basename(args.solute),
                solutename, args.molecule, nsolvent))

    if args.makendx :
        filename = os.path.splitext(args.outtop)[0]+".ndx"
        with open(filename, "w") as f :

            f.write("[ System ]\n")
            for idx, atom in enumerate(solvated.atoms, 1) :
                f.write("%4d "%idx)
                if idx % 15 == 0 : f.write("\n")
            f.write("\n")

            f.write("[ AA ]\n")
            for idx, atom in enumerate(solute.residues[0].atoms, 1) :
                f.write("%4d "%idx)
                if idx % 15 == 0 : f.write("\n")
            f.write("\n")

            f.write("[ VS ]\n")
            idx0 = len(solute.residues[0].atoms)
            idx = idx0
            for res in solute.residues[1:] :
                for atom in res.atoms :
                    idx += 1
                    f.write("%4d "%idx)
                    if (idx-idx0) % 15 == 0 : f.write("\n")
            f.write("\n")

            f.write("[ VSCG ]\n")
            idx0 = len(solute.residues[0].atoms)
            idx  = idx0
            for cgres in solvated.residues[1:] :
                for atom in cgres.atoms :
                    idx += 1
                    f.write("%4d "%idx)
                    if (idx-idx0) % 15 == 0 :
                        f.write("\n")
            f.write("\n")

            f.write("[ CG ]\n")
            idx0 = len(solute.atoms)
            idx  = idx0
            for cgres in solvated.residues[len(solute.residues):] :
                for atom in cgres.atoms :
                    idx += 1
                    f.write("%4d "%idx)
                    if (idx-idx0) % 15 == 0 :
                        f.write("\n")
            f.write("\n")
