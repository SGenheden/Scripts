# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to replace some POPI residues in a membrane patch with IPC lipids.
Assumes Slipid force field.

Used in membrane engineering project

Examples
--------
popimem2ipcmem.py -f ../32popi_22dopc_10erg/em_whole.gro -t ipc.pdb -z ipc.zmat -i ../slipids.ff/ipc.itp -n 10  > replaced.gro
"""
import sys
import argparse
import os

import numpy as np
import numpy.random as random

from sgenlib import pdb
from sgenlib import geo
from sgenlib import gmx

def convert_single(popi_xyz,ipc_xyz,ipc_zmat) :
    """
    Convert a single POPI lipid to a IPC lipid
    """
    # Will copy these atoms straight off
    atoms_to_copy = "C11 C12 C13 C14 C15 C16 H1 H2 H3 H4 H5 H6 O2 O3 O4 O5 O6 HO2 HO3 HO4 HO5 HO6 O12 P O14 O13 O11 C1 HA HB C2 HS".split()
    # Will overlay these atoms with POPI atoms with other names
    atoms_to_overlay = "NF  C1F C2F C3F C4F C5F C6F C7F C8F C9F C10F C11F C12F C13F C14F C15F C16F C17F C18F C3S C4S C5S C6S C7S C8S C9S C10S C11S C12S C13S C14S C15S C16S C17S C18S C19S C20S".strip().split()
    atoms_in_popi =    "O21 C21 C22 C23 C24 C25 C26 C27 C28 C29 C210 C211 C212 C213 C214 C215 C216 C217 C218 C3  O31 C31 C32 C33 C34 C35 C36  C37  C38  C39  C310 C311 C312 C313 C314 C315 C316".strip().split()
    # Create the correspondance list
    atoms_corr = {}
    for ipc_atom,popi_atom in zip(atoms_to_overlay,atoms_in_popi) :
        atoms_corr[ipc_atom] = popi_atom

    # Compute the value of the z-matrix from the pdb-file
    # This will be used to build atoms not taken from POPI
    zval = np.zeros([len(ipc_zmat),3])
    for i,atoms in enumerate(ipc_zmat) :
        if i < 3 : continue
        zval[i,0] = np.sqrt(np.sum((ipc_xyz[atoms[0]]-ipc_xyz[atoms[1]])**2))
        zval[i,1] = geo.angle_atms(ipc_xyz[atoms[0]],ipc_xyz[atoms[1]],ipc_xyz[atoms[2]])*180.0/np.pi
        zval[i,2] = geo.dihedral_protoms(ipc_xyz[atoms[0]],ipc_xyz[atoms[1]],ipc_xyz[atoms[2]],ipc_xyz[atoms[3]])*180.0/np.pi

    # First, copy the atoms of the head group
    xyz = {}
    for atom in atoms_to_copy :
        xyz[atom] = np.array(popi_xyz[atom],copy=True)

    # Second, overlay atoms in the sphingo/ester group and carbons in the tails
    for atom in atoms_to_overlay :
        xyz[atom] = np.array(popi_xyz[atoms_corr[atom]],copy=True)

    # Third, build the carbons in the long IPC chain from the z-matrix of the previous carbon
    # this prevents most bad geometries
    otheratoms = ["C%dF"%i for i in range(19,27)]
    for i,atoms1 in enumerate(ipc_zmat) :
        if atoms1[0] not in otheratoms : continue
        n = int(atoms1[0][1:3])
        for j,atoms2 in enumerate(ipc_zmat) :
              if atoms2[0] != "C%dF"%(n-1) : continue
              zval[i,0] = np.sqrt(np.sum((xyz[atoms2[0]]-xyz[atoms2[1]])**2))
              zval[i,1] = geo.angle_atms(xyz[atoms2[0]],xyz[atoms2[1]],xyz[atoms2[2]])*180.0/np.pi
              zval[i,2] = geo.dihedral_protoms(xyz[atoms2[0]],xyz[atoms2[1]],xyz[atoms2[2]],xyz[atoms2[3]])*180.0/np.pi
              xyz[atoms1[0]] = geo.build_xyz(xyz[atoms1[1]],xyz[atoms1[2]],xyz[atoms1[3]],zval[i,0],zval[i,1],zval[i,2])

    # Fourth, build the rest of the atoms from the pre-defined z-matrix and from values in the template PDB-file
    for i,atoms in enumerate(ipc_zmat) :
        if i < 3 or atoms[0] in atoms_to_copy or atoms[0] in atoms_to_overlay or atoms[0] in otheratoms : continue
        xyz[atoms[0]] = geo.build_xyz(xyz[atoms[1]],xyz[atoms[2]],xyz[atoms[3]],zval[i,0],zval[i,1],zval[i,2])

    return xyz

def replace(atoms,ipc_hash,ipc_zmat,ipc_atomlist,zdisp,natom,nres)  :
    """
    Replace a POPI molecule with an IPC and write it out in GRO-format
    """
    # Hash the POPI atoms
    popi_hash = {}
    for atom in atoms : popi_hash[atom.name.strip()] = atom.xyz

    # Use the convert algorithm
    ipc_xyz = convert_single(popi_hash, ipc_hash, ipc_zmat)

    # Check for bad geometry, i.e. a too bend tail
    angle = geo.angle_atms(ipc_xyz["C1F"],ipc_xyz["C13F"],ipc_xyz["C26F"])*180.0/np.pi
    sys.stderr.write("Angle for residue %d is %.3f \n"%(nres+1,angle))
    #if angle < 110 : raise Exception("Bad geometry")

    # Print out the IPC in GRO-format
    nres = nres + 1
    for i,atom in enumerate(ipc_atomlist) :
        natom = natom + 1
        print "%5d%5s%5s%5d%8.3f%8.3f%8.3f"%(nres,"IPC",atom,natom,ipc_xyz[atom][0]/10.0,ipc_xyz[atom][1]/10.0,ipc_xyz[atom][2]/10.0+zdisp)
    return natom,nres

if __name__ == '__main__':

    random.seed(1124142)

    # We will displace the upper part of the membrane this much to accomodate the large IPC chains
    Z_CONST = 1.5
    
    parser = argparse.ArgumentParser(description="Program to convert a POPI membrane to a membrane with IPC")
    parser.add_argument('-f','--file',help="the input gro-file")
    parser.add_argument('-t','--template',help="a template IPC pdb-file")
    parser.add_argument('-z','--zmat',help="the IPC z-matrix")
    parser.add_argument('-i','--itp',help="the IPC itp-file")
    parser.add_argument('-n','--nreplace',type=int,help="how many POPI molecules to replace")
    args = parser.parse_args()

    # Read in a gro-file
    popimem = pdb.PDBFile()
    popimem.read(args.file,gro=True)
    # Read in a template IPC-file
    ipc_pdb = pdb.PDBFile(filename=args.template)
    # Read in the z-matrix of the IPC molecule
    ipc_zmat  = [line.strip().split() for line in open(args.zmat,"r").readlines()]
    # Read in the atom order of the IPC itp-file
    ipctop = gmx.TopFile(args.itp)
    ipc_atomlist = [atom.name for atom in ipctop.moleculetypes[0].atoms]
    nreplace = args.nreplace

    # Check where the middle of the membrane is
    midz = popimem.xyz[:,2].mean()

    # Hash the IPC atoms
    ipc_hash = {}
    for atom in ipc_pdb.atoms : ipc_hash[atom.name.strip()] = atom.xyz

    # Find POPI residues in lower and upper part of the membrane
    popires_low = []
    popires_upp = []
    for residue in popimem.residues :
        if residue.resname[:3] == "POP" :
            if residue.collect("centerofmass")[2] > midz :
                popires_upp.append(residue)
            else :
                popires_low.append(residue)
    npopileft = len(popires_upp)+len(popires_low)-nreplace

    # Star of a GRO-file
    print "POPI replaced with IPC"
    print "%5d"%(len(popimem.atoms)-len(popires_low[0].atoms)*nreplace+len(ipc_zmat)*nreplace)

    natom = 0
    nres = 0

    # Replace nreplace/2 number of POPI molecules in the lower part of the membrane
    taken = [False]*len(popimem.residues)
    ntaken = 0
    while ntaken < nreplace/2 :
        # Find a residue that is not taken
        if npopileft > 0 :
            takei = random.randint(0,len(popires_low))
        else :
            takei = ntaken
        while taken[popires_low[takei].idx] : takei = random.randint(0,len(popires_low))
        # Try to replace the POPI molecule, with exceptions for bad geometries
        try :
            natom,nres=replace(popires_low[takei].atoms,ipc_hash,ipc_zmat,ipc_atomlist,0.0,natom,nres)
            taken[popires_low[takei].idx] = True
            ntaken = ntaken + 1
        except :
            pass
    sys.stderr.write("Ntaken = %d\n"%ntaken)

    # Replace nreplace/2 number of POPI molecules in the upper part of the membrane
    while ntaken < nreplace :
        if npopileft > 0 :
            takei = random.randint(0,len(popires_upp))
        else :
            takei = ntaken-len(popires_low)
        while taken[popires_upp[takei].idx] : takei = random.randint(0,len(popires_upp))
        try :
            natom,nres=replace(popires_upp[takei].atoms,ipc_hash,ipc_zmat,ipc_atomlist,Z_CONST,natom,nres)
            taken[popires_upp[takei].idx] = True
            ntaken = ntaken + 1
        except :
            pass
    sys.stderr.write("Ntaken = %d\n"%ntaken)

    # Write out all other residues in the membrane
    for i,residue in enumerate(popimem.residues) :
        if taken[i] : continue
        zdisp = 0.0
        if residue.collect("centerofmass")[2] > midz : zdisp = Z_CONST
        for atom in residue.atoms :
            natom = natom + 1
            print "%5d%5s%5s%5d%8.3f%8.3f%8.3f"%(i+nreplace,atom.resname,atom.name,natom,atom.x/10.0,atom.y/10.0,atom.z/10.0+zdisp)
    print "%8.3f%8.3f%8.3f"%(popimem.box[0]/10.0,popimem.box[1]/10.0,popimem.box[2]/10.0+Z_CONST)
