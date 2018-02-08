# Author: Samuel Genheden samuel.genheden@gmail.com
"""
Program to analyze hydrogen bonds
"""

import argparse
import cPickle

#import recsql
import numpy as np

"""
Head: O13, O14, O2-O6
Backbone: O22, O32, OF, O2F, O4S, O1S, NF
Sterol: OH2
Ergosterol: O3
"""

class Counts(object):
    def __init__(self, kind):
        self.heads = 0
        self.backbones = 0
        self.sterols = 0
        self.waters = 0
        self.ethanols = 0
        self.total = 0
        self.doubles = 0
        self.nsnap = 0
        self.nmol = None
        self.kind = kind

    def update(self, curr):
        self.nsnap += 1
        self.nmol = len(curr)
        for lst in curr:
            rlst = self._reduce(lst)
            self._classify(rlst)

    def _classify(self, row):
        if row is None : return

        self.total += 1
        res = getattr(row,self.kind+"_resnm")
        atom = getattr(row,self.kind+"_atom")
        #if self.kind == "donor" : print atom+res,
        if res == "ERG" :
            self.sterols += 1
        elif res == "SOL" :
            self.waters += 1
        elif res == "eth" :
            self.ethanols += 1
        elif atom in ["O13", "O14", "O2", "O3", "O4", "O5", "O6", "HO2", "HO3", "HO4", "HO5", "HO6"]:
            self.heads += 1
        elif atom in ["O22", "O32", "OF", "O2F", "O4S", "O1S", "NF", "HO2F", "HO4S", "HO3S", "HNF"]:
            self.backbones += 1

    def _reduce(self, lst):

        if len(lst) == 0 :
            return None
        elif len(lst) == 1 :
            return lst[0]
        else :
            self.doubles += 1
            idx = np.argmin([row.distance for row in lst])
            return lst[idx]

    def __str__(self):
        fac = 1.0  / float(self.nsnap*self.nmol)
        return "\n".join([
                "\tHead:\t%.4f"%(self.heads*fac),
                "\tBackbone:\t%.4f"%(self.backbones*fac),
                "\tSterol:\t%.4f"%(self.sterols*fac),
                "\tWater:\t%.4f"%(self.waters*fac),
                "\tEthanol:\t%.4f"%(self.ethanols*fac),
                "\tTotal:\t%.4f"%(self.total*fac),
                "\tDouble:\t%.4f"%(self.doubles*fac)
            ])

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description="Script analyse hbond db")
    argparser.add_argument('pickle', help="the pickle file")
    argparser.add_argument('--hac',type=int,help="the number of hac", default=3)
    argparser.add_argument('--ac',type=int,help="the number of ac", default=5)
    argparser.add_argument('--dt',type=float,help="the number of ps per snapshot", default=50.0)
    args = argparser.parse_args()

    nhac = args.hac
    nac = args.ac
    offset_hac = -5317 # We assume that the first hac is residue 5318
    offset_ac  = -5320

    p = cPickle.load(open(args.pickle))
    #q = recsql.SQLarray('tbl',p)
    donors_hac = Counts("acceptor")
    acceptors_hac = Counts("donor")
    if args.ac > 0:
        acceptors_ac = Counts("donor")

    prevtime = None
    curr_acceptors_hac = [[] for i in range(nhac)]
    curr_donors_hac = [[] for i in range(nhac)]
    if args.ac > 0:
        curr_acceptors_ac = [[] for i in range(nac)]
    for row in p:
        if prevtime is None or prevtime != row.time :
            donors_hac.update(curr_donors_hac)
            acceptors_hac.update(curr_acceptors_hac)
            curr_acceptors_hac = [[] for i in range(nhac)]
            curr_donors_hac = [[] for i in range(nhac)]
            if args.ac > 0:
                curr_acceptors_ac = [[] for i in range(nac)]
                acceptors_ac.update(curr_acceptors_ac)

        if row.donor_resnm == "AAC":
            curr_donors_hac[row.donor_resid+offset_hac].append(row)
        elif row.acceptor_resnm == "AAC":
            curr_acceptors_hac[row.acceptor_resid+offset_hac].append(row)
        elif row.acceptor_resnm == "ace" and args.ac > 0:
            curr_acceptors_ac[row.acceptor_resid+offset_ac].append(row)
        prevtime = row.time
        #if acceptors_hac.nsnap > 20 : break

    print "Donors (HAc):\n",donors_hac
    print "Acceptors (HAc):\n",acceptors_hac
    if args.ac > 0:
        print ""
        print "Acceptors (Ac-):\n",acceptors_ac
