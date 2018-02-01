# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to parse RESP charges and make Gromacs residue template file (.rtp)

Atoms in the PDB file need to be in the same order as in the charge file

The atom types file need to have an atomtype definition on each line
    NAME1 TYPE1
    NAME2 TYPE2
    ...

Used in membrane engineering project

Examples
--------
parse_optq.py -f model0_1.pdb -q qout -o model0.rtp -t atypes.txt
    Make an rtp file based on model0_1 and qout
"""
import argparse

import parmed

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description="Script to parse optimal charges")
    argparser.add_argument('-f','--file',help="the PDB file")
    argparser.add_argument('-q','--qout',help="the output charges",default="qout")
    argparser.add_argument('-o','--out',help="the output RTP file")
    argparser.add_argument('-t','--types',help="a file with atom types")
    args = argparser.parse_args()

    struct = parmed.load_file(args.file)

    qline = ""
    with open(args.qout, "r") as f :
        line = f.readline()
        while line :
            qline += line.strip() + " "
            line = f.readline()
    charges = map(float,qline.strip().split())

    for atom, charge in zip(struct.atoms, charges) :
        print "%4s%10.6f"%(atom.name, charge)

    if args.out is not None :
        atype = {}
        with open(args.types, "r") as f :
            for line in f.readlines() :
                a, t = line.strip().split()
                atype[a] = t

        with open(args.out, "w") as f :
            f.write("[ bondedtypes ]\n")
            f.write("1       5          9        2        1           3      1     0\n\n")
            f.write("[ UNK ]\n\n")
            f.write("[ atoms ]\n")
            for i, (atom, charge) in enumerate(zip(struct.atoms, charges)) :
                f.write("%5s %6s %10.6f %3d\n"%(atom.name,
                        atype[atom.name], charge, i))
            f.write("\n[ bonds ]\n")
            for bond in struct.bonds :
                f.write("%5s %5s\n"%(bond.atom1.name, bond.atom2.name))
            f.write("\n")
