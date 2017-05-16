
import argparse

import numpy as np

def _parse_martini(filename, dest) :

    atomtypes = []
    with open(filename, "r") as f :
        line = f.readline()
        while line  and line.find("[ atomtypes ]") == -1 :
            line = f.readline()

        line = f.readline()
        dest.write("; MARTINI atom types, including virtual types\n")
        while line and line[0] != "[" :
            if not (line[0] == ";" or len(line.strip().split()) == 0) :
                cols = line.strip().split()
                dest.write("%-9s %-8s %10.5f  %10.6f  A %13.6e %13.6e\n"%(
                    cols[0], cols[0], float(cols[1]), 0.0, 0.0, 0.0))
                dest.write("%-9s %-8s %10.5f  %10.6f  V %13.6e %13.6e\n"%(
                    "v"+cols[0], "v"+cols[0], 0.0, 0.0, 0.0, 0.0))
                atomtypes.append(cols[0])
            line = f.readline()

        dest.write("\n[ nonbond_params ]\n")
        dest.write(";%-8s %-8s   %10s  %10s \n"%(
            "itype", "jtype", "sigma", "epsilon"))
        dest.write("; Martini mixing rules\n")
        line = f.readline()
        while line and line[0] != "[" :
            if not (line[0] == ";" or len(line.strip().split()) == 0) :
                cols = line.strip().split()
                B = float(cols[3])
                A = float(cols[4])
                if B != 0.0 :
                    try :
                        sig = np.round((A/B)**(1.0/6.0), 2)
                        eps = np.round(B**2/(4*A), 2)
                    except :
                        raise Exception("Could not parse line: %s"%line)
                else :
                    sig = 0.0
                    eps = 0.0
                dest.write("%-9s %-8s 1 %10.5f  %10.5f\n"%(
                    cols[0], cols[1], sig, eps))
                dest.write("%-9s %-8s 1 %10.5f  %10.5f\n"%(
                    cols[0], "v"+cols[1], sig, eps))
                # V-V should be zero
                #dest.write("%-9s %-8s 1 %10.5f  %10.5f\n"%(
                #    "v"+cols[0], "v"+cols[1], sig, eps))
                dest.write("%-9s %-8s 1 %10.5f  %10.5f\n"%(
                    "v"+cols[0], cols[1], sig, eps))
            line = f.readline()
    return atomtypes

def _parse_gaff(filename, dest) :

    dest.write("[ atomtypes ]\n")
    dest.write(";%-8s %-8s %10s  %10s    %13s %13s\n"%(
        "type", "type", "mass", "charge", "sigma", "epsilon"))
    dest.write("; GAFF atom types\n")

    atomtypes = []
    with open(filename, "r") as f :
        line = f.readline()
        while line and line.find("MOD4      RE") == -1 :
            line = f.readline()

        line = f.readline()
        while line :
            col = line.strip().split()
            if len(col) < 3 : break
            atype = col[0]
            rmin = float(col[1])*0.1 # 0.1 to Convert to nm
            eps = float(col[2]) * 4.184 # To convert to kJ/mol
            dest.write("%-9s %-8s %10.5f  %10.6f  A %13.6e %13.6e ; %6.4f %6.4f\n"%(
                atype, atype, 0.0, 0.0, rmin * 2**(-1.0/6.0) * 2, eps ,
                rmin, eps))
            atomtypes.append(atype)
            line = f.readline()
    return atomtypes

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description="Script to make nonbonded parameters for hybrid MARTINI simulations")
    argparser.add_argument('-gaff','--gaff',help="the Amber gaff parameters",default="gaff.dat")
    argparser.add_argument('-martini','--martini',help="the martini parameters", default="martini_v2.2.itp")
    argparser.add_argument('-o', '--out', help="the output filename", default="ffnonbonded.itp")
    args = argparser.parse_args()

    with open(args.out, "w") as dest :
        dest.write("\n;Non-bonded parameters for hybrid GAFF/MARTINI simulations\n")
        dest.write("\n\n[ defaults ]\n")
        dest.write("; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ\n")
        dest.write("1               2               yes             0.5     0.8333\n\n")

        gafftypes = _parse_gaff(args.gaff, dest)
        martinitypes = _parse_martini(args.martini, dest)
        dest.write("\n\n;Hybrid LJ, only repulsive potential\n")
        for gtype in gafftypes :
            for mtype in martinitypes :
                dest.write("%-9s %-8s 1 %10.5e  %10.5e\n"%(
                    gtype, mtype, -0.1, 0.25e+05))
