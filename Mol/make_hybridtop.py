# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to create a hybrid MARTINI topology

The AA structure should have the same order as the structure
created by automartini. This is important since antechamber reshuffles
the atoms.

This script requires the branch write_virtualn of my version of parmed
available at https://github.com/SGenheden/ParmEd/

Examples:
    gmx_hybridtop.py -p AA=aa_gmx.top CG=cg_gmx.top -s aa.pdb -m map.dat
"""

import argparse

import parmed

def _load_mapping(filename, struct) :
    mapping = {}
    with open(filename, "r") as f :
         line = f.readline()
         while line:
            atom, atomlist = line.split(" : ")
            mapping[int(atom)] = [struct.atoms[int(i)-1].name
                                    for i in atomlist.strip().split()]
            line = f.readline()
    return mapping

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Program create a hybrid MARTINI topology")
    parser.add_argument('-p','--topol',nargs=2, help="the topology files")
    parser.add_argument('-s','--struct', help="the AA structure file")
    parser.add_argument('-m','--map', help="the CG to AA map")
    parser.add_argument('-o','--out',help="the output top-file",default="hybrid.itp")
    args = parser.parse_args()

    intop = {}
    param = {"AA":True, "CG":False}
    for t in args.topol :
        key, nam = t.split("=")
        intop[key] = parmed.gromacs.GromacsTopologyFile(nam, parametrize=param[key])

    # The CG atom will have a mass of zero as is required for virtual
    # sites and the atom type will be special as well, prepended with a v
    for atom in intop["CG"].atoms :
        atom.mass = 0.0
        atom.type = "v"+atom.type

    # Remove connectivity of CG topology as these beads
    # will be created as virtual sites
    intop["CG"].bonds = []
    intop["CG"].angles = []
    intop["CG"].dihedrals = []

    # Load an AA structures that was created from automartini,
    # this is necessary since antechamber has reshuffled the atoms
    aastruct = parmed.load_file(args.struct)
    mapping = _load_mapping(args.map, aastruct)

    # This combines the AA and CG topologies, neat right!
    hybridtop = intop["AA"]+intop["CG"]

    # This adds the virtual_sitesn directives to the topology file
    naa = len(intop["AA"].atoms)
    for atom, atomlist in mapping.iteritems() :
        obj_list = []
        for aa_atom in hybridtop.atoms[:naa] :
            if aa_atom.name in atomlist :
                obj_list.append(aa_atom)
        siteatom = hybridtop.atoms[naa + atom - 1]
        v = parmed.VirtualSiteN(siteatom, obj_list)
        hybridtop.virtual_sitesn.append(v)
    hybridtop.write(args.out, combine="all", itp=True)
