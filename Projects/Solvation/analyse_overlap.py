# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to analyse the overlap of solutes in two solvents
The solutes are taken from the Minnesota solvation database
"""

import argparse
import re

import dblib

def _parse_info(form) :
    """
    Parse weight, number of atoms and number of heavy atoms from molecular formula
    """
    w = 0
    n = 0
    nh = 0
    for part in re.findall("[A-Z]+[0-9]+",form):
        m = re.match("([A-Z]+)([0-9]+)",part)
        element = m.group(1)
        number = int(m.group(2))
        w += mass[element.capitalize()]*number
        n += number
        if element != "H" : nh += number
    return w,n,nh

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description="Script to analyse overlap of solutes in two solvents")
    argparser.add_argument('-db', '--db', help="the molecule database")
    argparser.add_argument('-solvent1', '--solvent1', help="the first solvent")
    argparser.add_argument('-solvent2', '--solvent2', help="the second solvent")
    args = argparser.parse_args()

    db = dblib.SolvDb(filename=args.db,type="abs",filehandle="^0")


    nfound = 0
    for entry1,entry2 in db.itersoluteoverlap(args.solvent1,args.solvent2) :
        w,n,nh = _parse_info(entry1.Formula)
        if (n-nh) > 2 and nh <= args.cutoff:
            print entry1.SoluteName+";"+entry1.Formula+";"+"%.3f"%w
            nfound += 1

        if args.plot :
            allw.append(w)
            alln.append(n)
            allnh.append(nh)

    print "____\n%d"%nfound

    if args.plot:
        for i,(data,name) in enumerate(zip([allw,alln,allnh],["allw","alln","allnh"]),1):
            f = plt.figure(i)
            h,e = np.histogram(np.asarray(data),bins=10)
            f.gca().plot((e[:-1]+e[1:])/2.0,h)
            f.savefig(name+".png",format="png")
