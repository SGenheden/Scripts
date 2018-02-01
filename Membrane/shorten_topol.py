# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to shorten the topology (itp-file) of a lipid such that
the acyl chains are shorter. Assumes Slipid force field.

Uses in membrane engineering project

Examples
--------
shorten_topol.py -f dopc.top -a 16 -o dopx.itp
    shorten DOPC to C16 chains
"""
import argparse
import os

import parmed

def _get_atom(parm, name) :
    """
    Return an Atom object with the specified name
    """
    for atom in parm.atoms  :
        if atom.name == name :
            return atom
    return None

def _remove_atoms(parm, last_idx, template, htemplate) :
    """
    Remove the terminal hydrogens, and make the terminal carbon a hydrogen,
    also re-distribute charges.

    Parameters
    ----------
    parm : parmed.Structure
        the topology
    last_idx : integer
        the index of the terminal carbon
    template : string
        the name prefix of the carbon atoms
    htemplate: string
        the name postfix of the hydrogen atoms
    """

    # Find all carbons preceeding the terminal carbon and the hydrogens
    # bonded to them
    carbons = []
    hydrogens = []
    for idx in range(last_idx - 2, last_idx + 1) :
        carbons.append(_get_atom(parm, "%s%d"%(template, idx)))
        hydrogens.append([_get_atom(parm, "H%d%s"%(idx, hstr))
                            for hstr in htemplate[:2]])

    # Start with removing the terminal hydrogens and its connectivity
    for hstr in htemplate :
        atom = _get_atom(parm, "H%d%s"%(last_idx, hstr))
        for bond in [bond for bond in parm.bonds if atom in bond] :
            parm.bonds.remove(bond)
        for angle in [angle for angle in parm.angles if atom in angle] :
            parm.angles.remove(angle)
        for ub in [ub for ub in parm.urey_bradleys if atom in ub] :
            parm.urey_bradleys.remove(ub)
        for dihedral in [dihedral for dihedral in parm.dihedrals if atom in dihedral] :
            parm.dihedrals.remove(dihedral)
        for pair in [pair for pair in parm.adjusts if atom in pair] :
            parm.adjusts.remove(pair)
        parm.atoms.remove(atom)

    # Change the terminal carbon to a hydrogen
    carbons[2].name = "H%d%s"%(last_idx-1,htemplate[2])
    carbons[2].type = "HAL3"
    carbons[2].charge = 0.016

    # Change the hydrogen and carbon names and type of the next to the terminal
    for h in hydrogens[1] :
        h.charge = 0.016
        h.type = "HAL3"
    carbons[1].charge = -0.081
    carbons[1].type = "CTL3"

    # Redistribute charges on hydrogen and carbons two positions from the terminal
    for h in hydrogens[0] :
        h.charge = -0.007
    carbons[0].charge = 0.047


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

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Program shorten the topology of lipid acyl chains")
    parser.add_argument('-f','--file',help="the input top-file")
    parser.add_argument('-o','--out',help="the output top-file",default="extended.top")
    parser.add_argument('-a','--atom',type=int,help="the new end atom")
    parser.add_argument('-c','--chains',nargs="+",help="the chains to extend",default=["sn1", "sn2"])
    parser.add_argument('-r','--resname',help="the new residue name")
    args = parser.parse_args()

    parm = parmed.load_file(args.file, parametrize=False)

    parm.residues[0].name = args.resname

    if "sn2" in args.chains :
        last_idx = _find_chain_end(parm, 20, "C2")
        for idx in range(last_idx, args.atom, -1) :
            _remove_atoms(parm, idx, "C2", ["R", "S", "T"])

    if "sn1" in args.chains :
        last_idx = _find_chain_end(parm, 20, "C3")
        for idx in range(last_idx, args.atom, -1) :
            _remove_atoms(parm, idx, "C3", ["X", "Y", "Z"])

    if os.path.splitext(args.out)[1] == ".itp" :
        parm.write(args.out, itp=True)
    else :
        parm.write(args.out)
