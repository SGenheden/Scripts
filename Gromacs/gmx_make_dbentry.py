# Author: Samuel Genheden, samuel.genheden@gmail.com

"""
Program to make a Gromacs residue template file
"""

import argparse

import parmed

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Program to make RTP entries")
    parser.add_argument('-f','--file',help="the input structure file")
    parser.add_argument('-t','--template',help="a template top-file")
    args = parser.parse_args()

    struct = parmed.load_file(arg.file)
    template = parmed.load_file(arg.template, parametrize=False)

    print """[ bondedtypes ]
1       5          9        2        1           3      1     0

[ UNK ]

[ atoms ]"""

    for i, atom in enumerate(struct) :
        try :
            templ_atom = template["@%s"%atom.name][0]
        except :
            templ_atom = None
        charge = templ_atom.charge if templ_atom is not None else 0.0
        atype = templ_atom.type if templ_atom is not None else "XXX"
        print "%5s %5s %10.6f %3d"%(atom.name, atype, charge, i)

    print "\n[ bonds ]"
    printed = []
    for atom1 in struct.atoms :
        for atom2 in atom1.bond_partners :
            if atom1.idx < atom2.idx :
                bond = "%5s %5s"%(atom1.name, atom2.name)
            else :
                bond = "%5s %5s"%(atom2.name, atom1.name)
            if bond not in printed :
                print bond
                printed.append(bond)
