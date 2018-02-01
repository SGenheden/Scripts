# Author: Samuel Genheden, samuel.genheden@gmail.com

"""
Program to setup LAMMPS restraints for APR calculations
"""

import argparse

import numpy as np

from sgenlib import geo
from sgenlib import pdb

colvarstr="""
colvarsTrajFrequency 2500
colvarsRestartFrequency 0

colvar {
    name d1
    width 1.0
    distance {
        group1 { atomNumbers 1 }
        group2 { atomNumbers %d }
    }
}

colvar {
    name a1
    width 1.0
    angle {
        group1 { atomNumbers 2 }
        group2 { atomNumbers 1 }
        group3 { atomNumbers %d }
    }
}

colvar {
    name a2
    width 1.0
    angle {
        group1 { atomNumbers 1 }
        group2 { atomNumbers %d }
        group3 { atomNumbers %d }
    }
}

harmonic {
    colvars d1
    forceConstant %.4f
    centers %.4f
}

harmonic {
    colvars a1
    forceConstant %.4f
    centers %.4f
}

harmonic {
    colvars a2
    forceConstant %.4f
    centers %.4f
}
"""

def atomindex(struct, ambmask) :

    atomnam = ambmask.split('@')[1]
    resname = ambmask.split('@')[0][1:]
    for i, atom in enumerate(struct.atoms) :
        if atom.name.strip() == atomnam and atom.resname.strip() == resname :
            return i
    return None

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description="Script to LAMMPS restraints for APR")
    argparser.add_argument('-i','--input',help="the apr.in file",default="apr.in")
    argparser.add_argument('-s','--structure',help="a structure file",default="align_z.pdb")
    argparser.add_argument('-p','--pull',type=int,help="the number of pull windows")
    argparser.add_argument('-d','--displ',type=float,nargs="+",help="the displacement vector")
    argparser.add_argument('-f','--force',type=float,nargs=2,help="the bond and angle force respectively",default=[5, 100])
    argparser.add_argument('-l','--label',help="a label for the output files",default="restraints")
    args = argparser.parse_args()

    struct = pdb.PDBFile(args.structure)

    h1 = ""
    h2 = ""
    h3 = ""
    g1 = ""
    g2 = ""
    with open(args.input, "r") as f :
        for line in f.readlines() :
            if line.strip().startswith("H1") :
                h1 = line.strip().split("=")[1].strip()
            elif line.strip().startswith("H2") :
                h2 = line.strip().split("=")[1].strip()
            elif line.strip().startswith("H3") :
                h3 = line.strip().split("=")[1].strip()
            elif line.strip().startswith("G1") :
                g1 = line.strip().split("=")[1].strip()
            elif line.strip().startswith("G2") :
                g2 = line.strip().split("=")[1].strip()
    print "# Atom selections:", h1, h2, h3, g1, g2

    h1 = atomindex(struct, h1)
    h2 = atomindex(struct, h2)
    h3 = atomindex(struct, h3)
    g1 = atomindex(struct, g1)
    g2 = atomindex(struct, g2)
    print "# Atom indices: ", h1, h2, h3, g1, g2

    x = struct.xyz
    eq1 = np.sqrt(np.sum((x[0, :]-x[h1,:])**2))
    eq2 = geo.angle_atms(x[1, :], x[0, :], x[h1,:])*180/np.pi
    eq3 = geo.dihedral(x[2, :], x[1, :], x[0, :], x[h1,:])*180/np.pi
    eq4 = geo.angle_atms(x[0, :], x[h1, :], x[h2,:])*180/np.pi
    eq5 = geo.dihedral(x[1, :], x[0, :], x[h1,:], x[h2,:])*180/np.pi
    eq6 = geo.dihedral(x[0, :], x[h1,:], x[h2,:], x[h3, :])*180/np.pi
    eq7 = np.sqrt(np.sum((x[0, :]-x[g1,:])**2))
    eq8 = geo.angle_atms(x[1, :], x[0, :], x[g1,:])*180/np.pi
    eq9 = geo.angle_atms(x[0, :], x[g1,:], x[g2,:])*180/np.pi

    with open("forcefield.%s_host"%args.label, 'w') as f :
        f.write("fix rest1 all restrain bond 1 "
            "%d ${bond_force0} ${bond_force0} %.4f\n"%(h1+1, eq1))
        f.write("fix rest2 all restrain angle 2 1 %d "
            "${angle_force0} ${angle_force0} %.4f\n"%(h1+1, eq2))
        f.write("fix rest3 all restrain dihedral 3 2 1 %d "
            "${angle_force0} ${angle_force0} %.4f\n"%(h1+1, -1*eq3+180))
        f.write("fix rest4 all restrain angle 1 %d %d "
            "${angle_force0} ${angle_force0} %.4f\n"%(h1+1, h2+1, eq4))
        f.write("fix rest5 all restrain dihedral 2 1 %d %d "
            "${angle_force0} ${angle_force0} %.4f\n"%(h1+1, h2+1, -1*eq5+180))
        f.write("fix rest6 all restrain dihedral 1 %d %d %d "
            "${angle_force0} ${angle_force0} %.4f\n"%(h1+1, h2+1, h3+1, -1*eq6+180))

    with open("forcefield.%s_guest"%args.label, 'w') as f :
        f.write("fix rest7 all restrain bond 1 %d "
            "${bond_force} ${bond_force} %.4f\n"%(g1+1, eq7))
        f.write("fix rest8 all restrain angle 2 1 %d "
            "${angle_force} ${angle_force} %.4f\n"%(g1+1, eq8))
        f.write("fix rest9 all restrain angle 1 %d %d "
            "${angle_force} ${angle_force} %.4f\n"%(g1+1, g2+1, eq9))

    # Conversions factors for colvars units
    angfac = (np.pi / 180.0) * (np.pi / 180.0)
    args.force = np.asarray(args.force) * 2
    if args.pull is not None :
        for displ in range(args.pull) :
            with open("colvars_p%d.in"%displ, "w") as f :
                f.write(colvarstr%(g1+1, g1+1, g1+1, g2+1, args.force[0], eq7+displ,
                    args.force[1]*angfac, eq8, args.force[1]*angfac, eq9))
    elif args.displ is not None  :
        for i, displ in enumerate(args.displ) :
            with open("%s_colvars_p%d.in"%(args.label,i), "w") as f :
                f.write(colvarstr%(g1+1, g1+1, g1+1, g2+1, args.force[0], eq7+displ,
                    args.force[1]*angfac, eq8, args.force[1]*angfac, eq9))
