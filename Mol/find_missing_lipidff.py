# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to find missing force field parameters and look them up in CHARMM36

Parses a log file from grompp and files in the toppar archive

Used in membrane engineering project

Examples
--------
find_missing_lipidff.py -l grompp.log -p model0.top  -f toppar/par_all36_cgenff.prm toppar/stream/lipid/toppar_all36_lipid_bacterial.str
"""
import argparse
import re

toKJ = 4.184

def _find_params(cons, ff, parsemulti=False) :
    """
    Finds the missing parameters

    Parameters
    ----------
    cons : list of strings
        the atoms in the bond, angle or dihedral
    ff : list of strings
        the CHARMM36 force field
    parsemulti : bool
        if we should look for more than one parameter, e.g. for dihedrals

    Returns
    -------
    list of strings
        the parsed parameter lines
    """
    result = []
    for con in cons :
        r1 = re.compile("^"+"\s+".join(con.split()))
        r2 = re.compile("^"+"\s+".join(con.split()[::-1]))
        found = False
        multi = []
        for line in ff :
            if r1.match(line) or r2.match(line):
                idx = len(con.split())
                ok = True
                try :
                    float(line.strip().split()[idx])
                except :
                    ok = False
                if ok :
                    found = True
                    print line
                    if not parsemulti :
                        result.append(line)
                        break
                    else :
                        multi.append(line)
        if not found :
            print "Not found: "+con
            result.append(None)
        else :
            if parsemulti :
                result.append(multi)
    return result

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description="Script to find missing lipid ff")
    argparser.add_argument('-l','--log',help="the grompp log file")
    argparser.add_argument('-p','--top',help="the top file")
    argparser.add_argument('-f','--ff',nargs="+",help="the force field files")
    args = argparser.parse_args()

    # Parse the grompp log file
    loglines = []
    with open(args.log, "r") as f :
        loglines = f.readlines()
    errlines = []
    for i, line in enumerate(loglines) :
        if line.startswith("  No default") :
            errline = loglines[i-1]
            one, two = errline.split("line ")
            errlines.append(int(two.split("]")[0]))

    # Parse the topology file
    toplines = []
    with open(args.top, "r") as f :
        toplines = f.readlines()

    atomtypes = []
    i = 0
    while toplines[i].find("[ atoms ]") == -1 :
        i +=1
    i += 1
    while toplines[i].find("[ bonds ]") == -1 :
        if (not toplines[i].startswith(";")) and len(toplines[i].strip()) > 5 :
            atomtypes.append(toplines[i].strip().split()[1])
        i +=1

    # Assemble all the missing bonds, angles and dihedrals
    bonds = set()
    angles = set()
    dihedrals = set()
    for err in errlines :
        cols = map(int,toplines[err-1].strip().split())
        atypes = [atomtypes[idx-1] for idx in cols[:-1]]
        if atypes[0] <= atypes[-1] :
            t = " ".join(atypes)
        else :
            t = " ".join(atypes[::-1])
        if len(cols) == 3 :
            bonds.add(t)
        elif len(cols) == 4 :
            angles.add(t)
        elif len(cols) == 5 :
            dihedrals.add(t)
    bonds = sorted(list(bonds), key=lambda atoms: atoms.split()[0])
    angles = sorted(list(angles), key=lambda atoms: atoms.split()[0])
    dihedrals = sorted(list(dihedrals), key=lambda atoms: atoms.split()[0])

    print "Missing bonds:"
    print "\n".join(bonds)
    print "\nMissing angles:"
    print "\n".join(angles)
    print "\nMissing dihedrals:"
    print "\n".join(dihedrals)

    # Read in all the parameter files
    fflines = []
    for filename in args.ff :
        with open(filename, "r") as f :
            fflines.extend([l.strip() for l in f.readlines()])

    print "\nFound bond parameters"
    bond_params = _find_params(bonds, fflines)
    print "\nFound angle parameters"
    angle_params = _find_params(angles, fflines)
    print "\nFound dihedral parameters"
    dihedral_params = _find_params(dihedrals, fflines, parsemulti=True)

    print "\nConverted bond parameters"
    for bondparm, bond in zip(bond_params, bonds) :
        if bondparm is None :
            print bond+" ???? "
            continue
        parm = bondparm.split("!")[0].strip().split()
        print "%-6s %-6s 1 %8.5f %10.2f"%(parm[0],parm[1],float(parm[3])*0.1,float(parm[2])*toKJ*100)

    print "\nConverted angle parameters"
    for angleparam, angle in zip(angle_params, angles) :
        if angleparam is None :
            print angle+" ???? "
            continue
        parm = angleparam.split("!")[0].strip().split()
        if len(parm) == 5 :
            print "%-6s %-6s %-6s 5 %7.2f %10.3f %8.5f %10.3f"%(parm[0],parm[1],
                parm[2],float(parm[4]),float(parm[3])*toKJ,0.0,0.0)
        else :
            print "%-6s %-6s %-6s 5 %7.2f %10.3f %8.5f %10.3f"%(parm[0],parm[1],
                parm[2],float(parm[4]),float(parm[3])*toKJ,float(parm[6])*0.1,
                float(parm[5])*toKJ*100)

    print "\nConverted dihedral parameters"
    for dihedparams, dihedral in zip(dihedral_params, dihedrals) :
        if dihedparams is None :
            print dihedral+" ???? "
            continue
        for dihedparam in dihedparams :
            parm = dihedparam.split("!")[0].strip().split()
            print "%-6s %-6s %-6s %-6s 9 %7.2f %10.3f %d"%(parm[0],parm[1], parm[2],
                parm[3],float(parm[6]),float(parm[4])*toKJ,int(parm[5]))
