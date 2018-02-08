# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to analyse chemical groups of solutes from the Minnesota solvation database

Examples:
analyse_chemicalgroups.py -db MNSol_alldata.txt -solvent hexanol -solutes hexanolwater.txt
"""

import argparse
import os
import subprocess
import tempfile

import dblib
from sgenlib import ambertools

babel_str = "babel -ipdb %s_leap.pdb -osdf %s_leap.sdf"
checkmol_str = "checkmol %s"

def _get_chemical_groups(filename):
    groups = []
    tmpfile,tmpname = tempfile.mkstemp()
    ret_code = subprocess.call(checkmol_str%filename, shell=True, stdout=tmpfile, stderr=tmpfile)
    # Catch some error codes
    if ret_code == 127:
        raise Exception("Unable to find checkmol executable, please make sure this is present in your PATH.")
    elif ret_code == 1:
        errmsg = "\n".join(line for line in open(tmpname).readlines())
        os.remove(tmpname)
        msg = "Check mol was not able to run successfully. Please check output. Error message was:\n%s"%(name,errmsg)
        raise Exception(msg)
    else:
        groups = [line.strip() for line in open(tmpname).readlines()]
        os.remove(tmpname)
    return groups

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description="Script to analyse chemical groups all solutes in a list")
    argparser.add_argument('-db', '--db', help="the molecule database")
    argparser.add_argument('-solvent', '--solvent', help="the solvent")
    argparser.add_argument('-solutes','--solutes',help="the list of solutes")
    args = argparser.parse_args()

    db = dblib.SolvDb(filename=args.db,type="abs",filehandle="^0")
    solutes = [s.strip() for s in open(args.solutes,'r').readlines()]

    groupset = set([])
    groupcount = {}
    sgroupcount = {}
    for entry in db.itersolutelist(args.solvent,solutes):
        try :
            with open(entry.FileHandle+".groups","r") as f :
                groups = [line.strip() for line in f.readlines()]
        except:
            ambertools.run_program("openbabel",babel_str%(entry.FileHandle,entry.FileHandle))
            groups = _get_chemical_groups(entry.FileHandle+"_leap.sdf")
            os.remove(entry.FileHandle+"_leap.sdf")
            with open(entry.FileHandle+".groups","w") as f :
                for group in groups :
                    f.write("%s\n"%group)

        if not groups :
            print "**",
            if entry.SoluteName.find("cyclo") == -1:
                groups.append("alkane")
            else:
                groups.append("cyclic hydrocarbons")
        sgroups = dblib.transform_groups(groups)
        print entry.SoluteName,":"," ".join(groups),":"," ".join(sgroups)
        groupset = groupset.union(set(groups))
        for group in groups :
            try :
                groupcount[group] += 1
            except :
                groupcount[group] = 1
        for group in sgroups:
            try:
                sgroupcount[group] += 1
            except:
                sgroupcount[group] = 1

    print ""
    for group in groupset:
        print group, groupcount[group]

    groupcount2 = {}
    for group in groupset :
        group2 = group.replace("primary","")
        group2 = group2.replace("secondary","")
        group2 = group2.replace("tertiary","")
        if group2.find("ether") > -1:
            group2 = "ether"
        elif group2.find("aliphatic amine") > -1:
            group2 = "amine"
        try :
            groupcount2[group2.strip()] += groupcount[group]
        except :
            groupcount2[group2.strip()] = groupcount[group]
    print ""
    for group,count in groupcount2.iteritems() :
        if count >= 5 : print group, count, sgroupcount[group]
