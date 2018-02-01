# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to replace a double bond in a lipid acyl chain with a cyclopropane group.
Assumes Slipid force field.

Used in membrane engineering project

Examples
--------
grow_cyclic.py -f single_popi.gro -r POPI -o single_popi_cyc.gro -t slipids.ff/popi.top -d popi_cyc.rtp
    Make a Gromacs residue template (rtp) file from a single lipid structure and topology
grow_cyclic.py -f 256lipids_md1_whole.gro -n 10 -r DOPC POPI -o m1_r4_cyc.gro --seed $RANDOM
    Replaces 10 random lipids in a membrane structure with cyclopropane lipids
"""
import argparse

import numpy as np
import numpy.random as random
import parmed
from parmed.geometry import _get_coords_from_atom_or_tuple

from sgenlib import geo

def _get_atom(parm, name) :
    """
    Return an Atom object with the specified name
    """
    for atom in parm.atoms  :
        if atom.name == name :
            return atom
    return None

def _get_atom_index(parm, name) :
    """
    Return the index of an Atom object with the specified name
    """
    for i, atom in enumerate(parm.atoms)  :
        if atom.name == name :
            return i
    return None

def _find_chain_end(parm, endidx, template) :
    """
    Find the index of the terminal carbon

    Parameters
    ----------
    parm : parmed.Structure
        the topology
    endidx : integer
        the very last index of that the terminal carbon can have
    template : string
        the name prefix of the carbon atoms

    Returns
    ---------
    integer
        the index
    """
    for idx in range(endidx-1, 2, -1) :
        a = _get_atom(parm, "%s%d"%(template, idx))
        if a is not None :
            return idx
    return -1

def _create_atoms(struct, residue, template, hnames) :
    """
    Adds the necessary atoms to a parmed.Structure object

    Parameters
    ----------
    struct : parmed.Structure
        the structure
    residue : string
        the residue name
    template : string
        the name prefix of the carbon atoms
    hnames: string
        the name postfox of the hydrogen atoms
    """
    try :
        h9 = struct.view["@H9%s&:%s"%(hnames[0],residue)].atoms[0]
    except :
        h9 = struct.view["@H9%s&:%s"%("1",residue)].atoms[0]
    c8 = struct.view["@%s8&:%s"%(template,residue)].atoms[0]
    c9 = struct.view["@%s9&:%s"%(template,residue)].atoms[0]
    c10 = struct.view["@%s10&:%s"%(template,residue)].atoms[0]
    term_i = _find_chain_end(struct.view[":%s"%residue], 30, template)
    terminal = struct.view["@%s%d&:%s"%(template,term_i,residue)].atoms[0]
    term_res_idx = _get_atom_index(terminal.residue, terminal.name)

    v1 = geo.vnorm(np.asarray(_get_coords_from_atom_or_tuple(c9)) -
                    np.asarray(_get_coords_from_atom_or_tuple(c10)))
    v2 = geo.vnorm(np.asarray(_get_coords_from_atom_or_tuple(c9)) -
                np.asarray(_get_coords_from_atom_or_tuple(h9)))
    v1_v2 = geo.vnorm(np.cross(v1,v2))
    xyz = 0.5*(np.asarray(_get_coords_from_atom_or_tuple(c9)) +
                    np.asarray(_get_coords_from_atom_or_tuple(c10))) + 1.3*v1_v2
    catom = parmed.Atom(list=struct.atoms, name="%s%d"%(template, term_i+1))
    catom.residue = c10.residue
    catom.xx = xyz[0]
    catom.xy = xyz[1]
    catom.xz = xyz[2]
    c10.residue.atoms.insert(term_res_idx+4, catom)
    struct.atoms.insert(terminal.idx+4, catom)
    struct.bonds.append(parmed.Bond(c9,catom))
    struct.bonds.append(parmed.Bond(c10,catom))

    xyz  = geo.build_xyz(np.asarray(_get_coords_from_atom_or_tuple(catom)),
                        np.asarray(_get_coords_from_atom_or_tuple(c9)),
                        np.asarray(_get_coords_from_atom_or_tuple(c8)),
                        0.95, 118, 150)
    hatom1 = parmed.Atom(list=struct.atoms, name="H%d%s"%(term_i+1, hnames[0]))
    hatom1.residue = catom.residue
    hatom1.xx = xyz[0]
    hatom1.xy = xyz[1]
    hatom1.xz = xyz[2]
    catom.residue.atoms.insert(term_res_idx+5, hatom1)
    struct.atoms.insert(terminal.idx+5, hatom1)
    struct.bonds.append(parmed.Bond(catom, hatom1))

    xyz  = geo.build_xyz(np.asarray(_get_coords_from_atom_or_tuple(catom)),
                        np.asarray(_get_coords_from_atom_or_tuple(c9)),
                        np.asarray(_get_coords_from_atom_or_tuple(c8)),
                        0.95, 118, 350)
    hatom2 = parmed.Atom(list=struct.atoms, name="H%d%s"%(term_i+1, hnames[1]))
    hatom2.residue = catom.residue
    hatom2.xx = xyz[0]
    hatom2.xy = xyz[1]
    hatom2.xz = xyz[2]
    catom.residue.atoms.insert(term_res_idx+6, hatom2)
    struct.atoms.insert(terminal.idx+6, hatom2)
    struct.bonds.append(parmed.Bond(catom, hatom2))

def _make_rtp(struct, template, mod_atoms, output) :
    """
    Make a Gromacs residue template

    Parameters
    ----------
    struct : parmed.Structure
        the structure
    template : parmed.Structure
        a structure to take charges and atom types from
    mod_atoms : list of strings
        the atoms that need to be modified with library parameters, not template
    output : string
        the filename of the RTP file
    """
    start_str ="""[ bondedtypes ]
1       5          9        2        1           3      1     0

[ UNK ]

[ atoms ]"""

    template = parmed.load_file(template, parametrize=False)
    libcharges = {}
    libcharges["9"] = -0.092305
    libcharges["10"] = -0.092305
    libcharges["8"] = -0.008999
    libcharges["11"] = -0.008999
    libcharges["19"] = -0.313038
    libcharges["H9"] = 0.083760
    libcharges["H10"] = 0.083760
    libcharges["H8"] = 0.026256
    libcharges["H11"] = 0.026256
    libcharges["H19"] = 0.121551

    libtypes = {}
    libtypes["9"] = "CG3RC1"
    libtypes["10"] = "CG3RC1"
    libtypes["19"] = "CG3C31"
    libtypes["H19"] = "HGA2"
    libtypes["H9"] = "HGA1"
    libtypes["H10"] = "HGA1"

    with open(output, "w") as f :
        f.write(start_str+"\n")
        for i, atom in enumerate(struct) :
            try :
                templ_atom = template["@%s"%atom.name][0]
            except :
                templ_atom = None
            charge = templ_atom.charge if templ_atom is not None else 0.0
            atype = templ_atom.type if templ_atom is not None else "XXX"

            if atom.name[0] == "C" :
                libname = atom.name[2:]
            elif atom.name[:-1] in ["H8","H11","H9", "H10","H19"] and atom.name[-1] in ["R", "S", "X", "Y","1"]:
                libname = atom.name[:-1]
            else :
                libname = "XXX"
            if libname in libcharges and atom.name in mod_atoms :
                charge = libcharges[libname]
                atype = libtypes[libname] if libname in libtypes else atype

            f.write("%5s %8s %10.6f %3d\n"%(atom.name, atype, charge, i))

        f.write("\n[ bonds ]\n")
        printed = []
        for atom1 in struct.atoms :
            for atom2 in atom1.bond_partners :
                if atom1.idx < atom2.idx :
                    bond = "%5s %5s"%(atom1.name, atom2.name)
                else :
                    bond = "%5s %5s"%(atom2.name, atom1.name)
                if bond not in printed :
                    f.write(bond+"\n")
                    printed.append(bond)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Program grow a cyclic propane structure")
    parser.add_argument('-f','--file',help="the input structure file")
    parser.add_argument('-o','--out',help="the output structure file")
    parser.add_argument('-r','--residue', nargs="+",help="the residue to adapt")
    parser.add_argument('-n','--nreplace',type=int,help="the number of residues to change")
    parser.add_argument('-t','--template',help="a template top-file")
    parser.add_argument('-d','--rtp',help="an rpt file")
    parser.add_argument('--seed', type=int, help="the random seed", default=1124142)
    args = parser.parse_args()

    random.seed(args.seed)

    struct = parmed.load_file(args.file)
    hnames = ["R", "S", "X", "Y"]
    mod_atoms = []

    if args.nreplace is None : # A single residue, will create a RTP file
        do_chain = _get_atom(struct, "H9S") is None
        if do_chain :
            print "Doing first chain"
            _create_atoms(struct, args.residue[0], "C2", hnames=hnames[:2])
            mod_atoms.extend(["C29","C210","C28","C211","C219","H9R","H10R","H91","H101","H8R","H8S","H11R","H11S","H19S","H19R"])

            do_chain = _get_atom(struct, "H9Y") is None
            if do_chain :
                print "Doing second chain"
                _create_atoms(struct, args.residue[0], "C3", hnames=hnames[2:])
                mod_atoms.extend(["C39","C310","C38","C311","C319","H9X","H10X","H8X","H8Y","H11X","H11Y","H19X","H19Y"])

        struct.save(args.out, overwrite=True)
        if args.rtp :
            _make_rtp(struct, args.template,  mod_atoms, args.rtp)

    else : # Will modify a membrane patch
        midz = struct.coordinates[:,2].mean()
        for resname in args.residue :
            residues = struct.view[":%s"%resname].residues
            res_low = []
            res_upp = []
            taken = {residue : False for residue in residues}
            # Find lipids in upper and lower leaflet
            for residue in residues :
                patom = None
                for atom in residue :
                    if atom.name == "P" :
                        patom = atom
                if patom.xz > midz :
                    res_upp.append(residue)
                else :
                    res_low.append(residue)
            print "Will replace %d %s residues"%(args.nreplace, resname)

            ntaken = 0
            while ntaken < int(np.ceil(args.nreplace * 0.5)) :
                takei = random.randint(0, len(res_low))
                while taken[res_low[takei]] :
                    takei = random.randint(0, len(res_low))
                if _get_atom(res_low[takei], "H9S") is None :
                    _create_atoms(struct, res_low[takei].number, "C2", hnames=hnames[:2])
                if _get_atom(res_low[takei], "H9Y") is None :
                    _create_atoms(struct, res_low[takei].number, "C3", hnames=hnames[2:])
                res_low[takei].name = res_low[takei].name[0]+"C"+res_low[takei].name[2:]
                taken[res_low[takei]] = True
                ntaken += 1

            while ntaken < args.nreplace :
                takei = random.randint(0, len(res_upp))
                while taken[res_upp[takei]] :
                    takei = random.randint(0, len(res_upp))
                if _get_atom(res_upp[takei], "H9S") is None :
                    _create_atoms(struct, res_upp[takei].number, "C2", hnames=hnames[:2])
                if _get_atom(res_upp[takei], "H9Y") is None :
                    _create_atoms(struct, res_upp[takei].number, "C3", hnames=hnames[2:])
                res_upp[takei].name = res_upp[takei].name[0]+"C"+res_upp[takei].name[2:]
                taken[res_upp[takei]] = True
                ntaken += 1

        struct.save(args.out, overwrite=True)
