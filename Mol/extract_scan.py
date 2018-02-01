# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to extract information from a Gaussian dihedral scan

It will write out energies and plot them to a .png file

I can create input for single-point calculations and PDB files for future use

Used in membrane engineering project

Examples
--------
extract_scan.py theta1_opt.log -c 6 5 1 2 -p "theta1_mp2=#MP2/cc-pVQZ" "theta1_cc=#CCSD(T)/cc-pVDZ" -s model_t1.pdb
    Extract from theta1_opt.log and create new output files for single-point optimisation with two energy functions
"""
import argparse
import os

import numpy as np
import matplotlib.pylab as plt
import parmed

from sgenlib import geo
from sgenlib import colors

AU2KJMOL=2625.5002

def _extract_energies(loglines, energy) :
    """
    Extract the energies from the Gaussian log file
    """
    i = len(loglines) - 8
    line = ""
    while not loglines[i].startswith(" N-N") :
        line = loglines[i].strip()+line
        i -= 1
    for data in line.split("\\") :
        if data.startswith(args.energy) :
            idx = len(args.energy) + 1
            energiesline = data[idx:]
            break
    return np.asarray(map(float,energiesline.split(",")))*AU2KJMOL

def _extract_coordinates(traj, coordinate) :
    """
    Extract the coordinates of all optimized steps
    """
    scanvalues = []
    for atoms in traj :
        if len(coordinate) == 4 :
            a1 = np.asarray(atoms[coordinate[0]-1][1:])
            a2 = np.asarray(atoms[coordinate[1]-1][1:])
            a3 = np.asarray(atoms[coordinate[2]-1][1:])
            a4 = np.asarray(atoms[coordinate[3]-1][1:])
            scanvalues.append(geo.dihedral(a1, a2, a3, a4)*180.0/np.pi)
    return scanvalues

def _setup_process(traj, label, job, nproc, charge, multiplicity) :
    """
    Create input for single point calculations
    """
    with open("%s.in"%label, "w") as f :
        for i, atoms in enumerate(traj, 1) :
            f.write("--Link1--\n")
            f.write("%%Chk=%s\n"%label)
            if nproc > 0 :
                f.write("%%NProcShared=%d\n"%nproc)
            jobline = job if i == 1 else job + " Guess=Read"
            f.write(jobline+"\n\n")
            f.write("%s coordinate set #%d\n\n"%(label,i))
            f.write("%s %s\n"%(charge, multiplicity))
            for atom in atoms :
                f.write("%5d%12.6f%12.6f%12.6f\n"%atom)
            f.write("\n\n")

def _write_xyz(loglines, filename) :
    """
    Write xyz file with all the optimized coordinates
    """
    traj = []
    i = 0
    while i < len(loglines) :
        while i < len(loglines) and loglines[i].find(" Optimization completed") == -1 :
            i += 1
        while i < len(loglines) and loglines[i].find("Standard orientation") == -1 :
            i += 1
        if i == len(loglines) :
            break
        i += 5
        atoms = []
        while loglines[i].strip()[0] != "-" :
            cols = loglines[i].strip().split()
            atoms.append((int(cols[1]),float(cols[3]),float(cols[4]),float(cols[5])))
            i += 1
        traj.append(atoms)
    atoms = traj.pop()
    traj.insert(0, atoms)

    with open(filename, "w") as f :
        for i, atoms in enumerate(traj, 1) :
            f.write("%d\n"%len(atoms))
            f.write("%d\n"%i)
            for atom in atoms :
                f.write("%5d%12.6f%12.6f%12.6f\n"%atom)

    return traj

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description="Script to extract information from Gaussian scan")
    argparser.add_argument('log',help="the Gaussian log file")
    argparser.add_argument('-s','--structure',help="a pdb structure")
    argparser.add_argument('-e','--energy',help="the energy flag", default="MP2")
    argparser.add_argument('-c','--coordinate',type=int,nargs="+",help="the coordinate that was scanned")
    argparser.add_argument('-p','--postprocess',nargs="+",help="setup postprocess jobs")
    args = argparser.parse_args()

    loglines = []
    with open(args.log, "r") as f :
        loglines = f.readlines()

    traj = _write_xyz(loglines, os.path.splitext(args.log)[0]+".xyz")

    energies = _extract_energies(loglines, args.energy)
    scanvalues = _extract_coordinates(traj, args.coordinate)
    print "%8s%12s"%("Scan value","Energy")
    for s, e in zip(scanvalues, energies) :
        print "%10.3f %12.3f"%(s, e)
    plt.plot(scanvalues, energies, "-*",color=colors.color(0))
    plt.xlabel("Scan value")
    plt.ylabel("Energy [kJ/mol]")
    plt.savefig(os.path.splitext(args.log)[0]+".png", dpi=360)

    if args.structure is not None :
        from zmat import mol
        graph = mol.GraphStructure(verbosity=0)
        graph.initialize(args.structure)
        graph.traverse()
        atoms = [z.split()[0] for z in graph.zmat]
        idx = np.asarray([atoms.index(atom.name) for atom in graph.structure])
        for i, atoms in enumerate(traj, 1) :
            filename = os.path.splitext(args.log)[0]+"_%d.pdb"%i
            graph.structure.coordinates = np.asarray([a[1:] for a in atoms])[idx]
            graph.structure.save(filename, overwrite=True)

    if args.postprocess is not None :
        nproc = 0
        for line in loglines :
            if line.find("NProcShared") > -1 :
                nproc = int(line.strip().split("=")[1])
                break
        charge = 0
        multiplicity = 1
        for line in loglines :
            if line.find("Charge = ") > -1 :
                charge = line.strip().split()[2]
                multiplicity = line.strip().split()[5]
                break
        for post in args.postprocess :
            lbl, job = post.split("=")
            _setup_process(traj, lbl, job, nproc, charge, multiplicity)
