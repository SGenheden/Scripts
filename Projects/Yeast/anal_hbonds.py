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

        return "".join([
                "\t%.4f"%(self.heads*fac),
                "\t%.4f"%(self.backbones*fac),
                "\t%.4f"%(self.sterols*fac),
                "\t%.4f"%(self.waters*fac),
                "\t %.4f"%(self.total*fac),
                "\t%.4f"%(self.doubles*fac)
            ])

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description="Script analyse hbond db")
    argparser.add_argument('--ethanols',type=int,help="the number of ethanols")
    argparser.add_argument('--dt',type=float,help="the number of ps per snapshot", default=10.0)
    args = argparser.parse_args()

    offset = -129 # We assume that the first ethanol is residue 129

    p = cPickle.load(open("mem_hbonds.pickle"))
    #q = recsql.SQLarray('tbl',p)
    donors = Counts("acceptor")
    acceptors = Counts("donor")

    prevtime = None
    curr_acceptors = [[] for i in range(args.ethanols)]
    curr_donors = [[] for i in range(args.ethanols)]
    for row in p:
        if prevtime is None or prevtime != row.time :
            donors.update(curr_donors)
            acceptors.update(curr_acceptors)
            curr_acceptors = [[] for i in  range(args.ethanols)]
            curr_donors = [[] for i in range(args.ethanols)]

        if row.donor_resnm in ["eth","but"]:
            curr_donors[row.donor_resid+offset].append(row)
        elif row.acceptor_resnm in ["eth","but"]:
            curr_acceptors[row.acceptor_resid+offset].append(row)
        prevtime = row.time
        #if row.time > 100 : break

    nsnap = row.time / args.dt
    print "Donors:\t",donors
    print "Acceptors:\t",acceptors
