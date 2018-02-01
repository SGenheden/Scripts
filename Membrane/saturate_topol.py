# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to convert a unsaturated acyl chain to a saturated chain in a
topology (itp-file). Works for Slipid force field.

Uses in membrane engineering project

Examples
--------
saturate_topol.py -f dopc.top -a C29 -o dopc_saturated.itp
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

def _new_angle(old, a1, a2) :
    """
    Create a new Angle object but one that replaces a1 with a2
    """
    if old.atom1 == a1 :
        new = parmed.Angle(a2, old.atom2, old.atom3)
    elif old.atom2 == a1 :
        new = parmed.Angle(old.atom1, a2, old.atom3)
    elif old.atom3 == a1 :
        new = parmed.Angle(old.atom1, old.atom2, a2)
    new.funct = 5
    return new

def _new_dihedral(old, a1, a2) :
    """
    Create a new Dihedral object but one that replaces a1 with a2
    """
    if old.atom1 == a1 :
        new = parmed.Dihedral(a2, old.atom2, old.atom3, old.atom4)
    elif old.atom2 == a1 :
        new = parmed.Dihedral(old.atom1, a2, old.atom3, old.atom4)
    elif old.atom3 == a1 :
        new = parmed.Dihedral(old.atom1, old.atom2, a2, old.atom4)
    elif old.atom4 == a1 :
        new = parmed.Dihedral(old.atom1, old.atom2, old.atom3, a2)
    new.funct = 9
    return new

def _new_atom(parm, carbon, hydrogen, name) :
    """
    Creates a new hydrogen atom and add the necessary topology

    Parameters
    ----------
    parm : parmed.Structure
        the topology
    carbon : parmed.Atom
        the carbon atom to which the hydrogen should be bonded
    hydrogen : parmed.Atom
        the other hydrogen, already in the topology
    name : string
        the name of the hydrogen atom

    Returns
    -------
    parmed.Atom
        the newly created atom
    """
    new = parmed.Atom(list=parm.atoms, name=name, type="HAL2", charge=0.0, mass=1.008000)
    new.residue = parm.residues[0]
    parm.residues[0].atoms.insert(hydrogen.idx+1, new)
    parm.atoms.insert(hydrogen.idx+1, new)
    # Add topology
    parm.bonds.append(parmed.Bond(carbon, new))
    for a in parm.adjusts :
        if a.atom1 == hydrogen :
            parm.adjusts.append(parmed.NonbondedException(new, a.atom2))
        elif a.atom2 == hydrogen :
            parm.adjusts.append(parmed.NonbondedException(a.atom1, new))
    for angle in hydrogen.angles :
        parm.angles.append(_new_angle(angle, hydrogen, new))
    parm.angles.append(parmed.Angle(hydrogen, carbon, new))
    parm.angles[-1].funct = 5
    for dihedral in hydrogen.dihedrals :
        parm.dihedrals.append(_new_dihedral(dihedral, hydrogen, new))
    return new

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Program saturate the topology of a lipid acyl chain")
    parser.add_argument('-f','--file',help="the input top-file")
    parser.add_argument('-o','--out',help="the output top-file",default="saturated.top")
    parser.add_argument('-a','--atom',help="the first atom of the double bond")
    args = parser.parse_args()

    parm = parmed.load_file(args.file, parametrize=False)

    template = args.atom[:2]
    idx = int(args.atom[2:])
    htemplate = ["R", "S"] if template[1] == "2" else ["X", "Y"]

    # Set the charges to zero for the surrounding carbon and hydrogen atoms
    # the atoms surrounding the double bond have charges, not the saturated chain
    _get_atom(parm, "%s%d"%(template, idx-1)).charge = 0.0
    _get_atom(parm, "H%d%s"%(idx-1, htemplate[0])).charge = 0.0
    _get_atom(parm, "H%d%s"%(idx-1, htemplate[1])).charge = 0.0
    _get_atom(parm, "%s%d"%(template, idx+2)).charge = 0.0
    _get_atom(parm, "H%d%s"%(idx+2, htemplate[0])).charge = 0.0
    _get_atom(parm, "H%d%s"%(idx+2, htemplate[1])).charge = 0.0

    # These are the atoms involved in the double bond
    a1 = _get_atom(parm, args.atom)
    h1 = _get_atom(parm, "H%d1"%idx)
    a2 = _get_atom(parm, "%s%d"%(template, idx+1))
    h2 = _get_atom(parm, "H%d1"%(idx+1))

    a1.charge = a2.charge = h1.charge = h2.charge = 0.0
    a1.type = a2.type = "CTL2"
    h1.type = h2.type = "HAL2"
    h1.name = "H%d%s"%(idx, htemplate[0])
    h2.name = "H%d%s"%(idx+1, htemplate[0])

    h1b = _new_atom(parm, a1, h1, "H%d%s"%(idx, htemplate[1]))
    h2b = _new_atom(parm, a2, h2, "H%d%s"%(idx+1, htemplate[1]))

    if os.path.splitext(args.out)[1] == ".itp" :
        parm.write(args.out, itp=True)
    else :
        parm.write(args.out)

    """print [a.name for a in _get_atom(parm, "H%d%s"%(idx-1, htemplate[0])).dihedral_partners]
    print [a.name for a in h1.dihedral_partners]
    print [a.name for a in h1b.dihedral_partners]
    print [a.name for a in h2.dihedral_partners]
    print [a.name for a in h2b.dihedral_partners]"""
