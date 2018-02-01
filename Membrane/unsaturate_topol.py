# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to make a polyunsaturated lipid, e.g. add an additional double bond in a
DOPC topology (itp-file). Works for Slipid force field.

Uses in membrane engineering project

Examples
--------
unsaturate_topol.py -f dopc.top -sn1 9 12 -sn2 9 12 -o dupc.itp
unsaturate_topol.py -f dopc.top -sn1 9 122 -sn2 9 12 15 -o llpc.itp
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

def _remove_atom(parm, hname) :
    """
    Remove a hydrogen and all connectivity associated with it

    Parameters
    ----------
    parm : parmed.Structure
        the topology
    hname : string
        the name of the hydrogen to remove
    """

    atom = _get_atom(parm, hname)
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

def _update_types(parm, startidx, ctemplate, hpostfix) :
    """
    Update the atom type to a "double bond"-type

    Parameters
    ----------
    parm : parmed.Structure
        the topology
    startidx : integer
        the index of the atom at the start of the double bond
    ctemplate : string
        the name template for the carbons
    hpostfix : string
        the final character of the hydrogen name
    """
    _get_atom(parm, "%s%d"%(ctemplate, startidx)).type = "CEL1"
    _get_atom(parm, "%s%d"%(ctemplate, startidx+1)).type = "CEL1"
    _get_atom(parm, "H%d%s"%(startidx, hpostfix)).type = "HEL1"
    _get_atom(parm, "H%d%s"%(startidx+1, hpostfix)).type = "HEL1"

def _update_charges(parm, startidx, ctemplate, hpostfix) :
    """
    Update the charges of a polyunsaturated chain

    Parameters
    ----------
    parm : parmed.Structure
        the topology
    startidx : list of integer
        the indices of the atom at the start of the double bonds
    ctemplate : string
        the name template for the carbons
    hpostfix : list of string
        the final characters of the hydrogen name
    """

    def _between_doubles(idx, startidx) :
        endidx = [i+1 for i in startidx]
        for end, start in zip(endidx[:-1],startidx[1:]) :
            if idx == end + 1 and idx == start -1 :
                return True
        return False

    for idx in range(startidx[0]-1,startidx[-1]+3) :
        catom = _get_atom(parm, "%s%d"%(ctemplate,idx))
        if catom.type == "CEL1" :
            catom.charge = -0.28
            _get_atom(parm, "H%d%s"%(idx, hpostfix[0])).charge = 0.14
        else :
            if _between_doubles(idx, startidx) :
                catom.charge = 0.28
                _get_atom(parm, "H%d%s"%(idx, hpostfix[0])).charge = 0.0
                _get_atom(parm, "H%d%s"%(idx, hpostfix[1])).charge = 0.0
            else :
                catom.charge = 0.12
                _get_atom(parm, "H%d%s"%(idx, hpostfix[0])).charge = 0.01
                _get_atom(parm, "H%d%s"%(idx, hpostfix[1])).charge = 0.01

    # Check if we have polyunsaturated chain ending close to the terminal
    # of the chain, there are special charges for this situation
    end_idx = _find_chain_end(parm, 30, ctemplate)
    if startidx[-1]+3 == end_idx :
        _get_atom(parm, "%s%d"%(ctemplate, end_idx)).charge = -0.27
        for postfix in hpostfix :
            _get_atom(parm, "H%d%s"%(end_idx, postfix)).charge = 0.05
        _get_atom(parm, "%s%d"%(ctemplate, end_idx-1)).charge = 0.26
        for postfix in hpostfix[:-1] :
            _get_atom(parm, "H%d%s"%(end_idx-1, postfix)).charge = 0.00


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Program unsaturate the topology of a lipid acyl chain")
    parser.add_argument('-f','--file',help="the input top-file")
    parser.add_argument('-o','--out',help="the output top-file",default="unsaturated.top")
    parser.add_argument('-sn1','--sn1',type=int,nargs="+",help="the double bonds in the sn1 chain")
    parser.add_argument('-sn2','--sn2',type=int,nargs="+",help="the double bonds in the sn2 chain")
    parser.add_argument('-r','--resname',help="the new residue name")
    args = parser.parse_args()

    parm = parmed.load_file(args.file, parametrize=False)

    parm.residues[0].name = args.resname

    if args.sn2 is not None :
        # First we will remove the hydrogen atoms and update the atom types
        for idx in args.sn2[1:] :
            _remove_atom(parm, "H%dS"%idx)
            _remove_atom(parm, "H%dS"%(idx+1))
            _update_types(parm, idx, "C2", "R")
        # Then we will update the charges
        _update_charges(parm, args.sn2, "C2", ["R", "S", "T"])


    if args.sn1 is not None :
        for idx in args.sn1[1:] :
            _remove_atom(parm, "H%dY"%idx)
            _remove_atom(parm, "H%dY"%(idx+1))
            _update_types(parm, idx, "C3", "X")
        _update_charges(parm, args.sn1, "C3", ["X", "Y", "Z"])

    if os.path.splitext(args.out)[1] == ".itp" :
        parm.write(args.out, itp=True)
    else :
        parm.write(args.out)
