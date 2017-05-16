# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to analyse chemical groups of solutes

Requires that the program checkmol is downloaded and installed

Examples:
chemical_groups.py -f mol1.sdf mol2.sdf
chemical_groups.py -f mol1.sdf mol2.sdf -o Groups/
"""

import argparse
import os
import subprocess
import tempfile

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

    argparser = argparse.ArgumentParser(description="Script to analyse chemical groups of solutes")
    argparser.add_argument('-f', '--files', nargs="+", help="the molecule files in sdf-format")
    argparser.add_argument('-o', '--out', help="the output folder")
    args = argparser.parse_args()

    for filename in args.files :
        groups = _get_chemical_groups(filename)

        dir = args.out if args.out is not None else os.path.dirname(filename)
        base = os.path.splitext(os.path.basename(filename))[0]
        filename = os.path.join(dir, base+".groups")
        with open(filename, "w") as f:
            for g in groups :
                f.write(g+"\n")
