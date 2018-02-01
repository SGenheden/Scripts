# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to extend the topology (itp-file) of a lipid such that
the acyl chains are longer. Assumes Slipid force field.

Used in membrane engineering project

Examples
--------
extend_topol.py -f dopc.top -a 20 -o dopx.itp
    extends DOPC to C20 chains
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

def _add_atom(parm, last_idx, template, htemplate) :
    """
    Change the terminal hydrogen to a carbon and add three hydrogen atoms,
    also add bonded information

    Parameters
    ----------
    parm : parmed.Structure
        the topology
    last_idx : integer
        the index of the terminal carbon
    template : string
        the name prefix of the carbon atoms
    htemplate: string
        the name postfox of the hydrogen atoms
    """

    # Find all carbons preceeding the terminal carbon and the hydrogens
    # bonded to them
    carbons = []
    hydrogens = []
    for idx in range(last_idx - 2, last_idx + 1) :
        carbons.append(_get_atom(parm, "%s%d"%(template, idx)))
        hydrogens.append([_get_atom(parm, "H%d%s"%(idx, hstr))
                            for hstr in htemplate[:2]])

    # Change the charges of the carbon and hydrogens next to the terminal carbon
    carbons[1].charge = 0.0
    hydrogens[1][0].charge = 0.0
    hydrogens[1][1].charge = 0.0

    # Change the type and charge of the terminal carbon and its hydrogens
    carbons[2].charge = 0.0470
    carbons[2].type = "CTL2"
    for h in hydrogens[2] :
        h.charge = -0.0070
        h.type = "HAL2"

    # Find the final hydrogen and change it to a carbon
    hatom = _get_atom(parm, "H%d%s"%(last_idx, htemplate[2]))
    hatom.type = "CTL3"
    hatom.charge = -0.0810
    hatom.name = "%s%d"%(template, last_idx+1)
    hatom.mass = 12.0110

    # Create 3 new hydrogen atoms and the topology information
    newatoms = []
    for i, hstr in enumerate(htemplate) :
        newatom = parmed.Atom(list=parm.atoms, name="H%d%s"%(last_idx+1, hstr), type="HAL3",
                                charge=0.0160, mass=1.008000)
        newatom.residue = parm.residues[0]
        parm.residues[0].atoms.insert(hatom.idx+1+i, newatom)
        parm.atoms.insert(hatom.idx+1+i, newatom)
        parm.bonds.append(parmed.Bond(hatom, newatom))
        parm.adjusts.append(parmed.NonbondedException(carbons[1], newatom))
        parm.adjusts.append(parmed.NonbondedException(hydrogens[2][0], newatom))
        parm.adjusts.append(parmed.NonbondedException(hydrogens[2][1], newatom))
        parm.angles.append(parmed.Angle(carbons[2], hatom, newatom))
        parm.dihedrals.append(parmed.Dihedral(carbons[1], carbons[2], hatom, newatom))
        parm.dihedrals.append(parmed.Dihedral(hydrogens[2][0], carbons[2], hatom, newatom))
        parm.dihedrals.append(parmed.Dihedral(hydrogens[2][1], carbons[2], hatom, newatom))
        newatoms.append(newatom)

    # Add angles between thre three added hydrogens
    parm.angles.append(parmed.Angle(newatoms[0], hatom, newatoms[1]))
    parm.angles.append(parmed.Angle(newatoms[0], hatom, newatoms[2]))
    parm.angles.append(parmed.Angle(newatoms[1], hatom, newatoms[2]))

    # Correct the angle and dihedral types
    for angle in parm.angles[-12:]:
        angle.funct = 5
    for dihedral in parm.dihedrals[-9:] :
        dihedral.funct = 9

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

    parser = argparse.ArgumentParser(description="Program extend the topology of lipid acyl chains")
    parser.add_argument('-f','--file',help="the input top-file")
    parser.add_argument('-o','--out',help="the output top-file",default="extended.top")
    parser.add_argument('-a','--atom',type=int,help="the new end atom")
    parser.add_argument('-c','--chains',nargs="+",help="the chains to extend",default=["sn1", "sn2"])
    parser.add_argument('-r','--resname',help="the new residue name")
    args = parser.parse_args()

    parm = parmed.load_file(args.file, parametrize=False)

    parm.residues[0].name = args.resname

    if "sn2" in args.chains :
        last_idx = _find_chain_end(parm, args.atom, "C2")
        for idx in range(last_idx, args.atom) :
            _add_atom(parm, idx, "C2", ["R", "S", "T"])

    if "sn1" in args.chains :
        last_idx = _find_chain_end(parm, args.atom, "C3")
        for idx in range(last_idx, args.atom) :
            _add_atom(parm, idx, "C3", ["X", "Y", "Z"])

    if os.path.splitext(args.out)[1] == ".itp" :
        parm.write(args.out, itp=True)
    else :
        parm.write(args.out)
