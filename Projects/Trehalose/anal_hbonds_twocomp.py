# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to analyse hydrogen bonds in a two-component system

It will classifies hydrogen bonds as trehalose-water (T-W),
trehalose-trehalose (T-T) or water-water (W-W). And for each
donor or acceptor atom it will write out the average number of hydrogen
bonds in each class.

Recognized donor atoms: H2O H3O H4O H6O HW1 HW2
Recognized acceptor atoms: O2 O3 O4 O6 OW

Examples:
    python anal_hbonds_twocomp.py --natoms 10262 -f r1_hbonds.pickle
    python anal_hbonds_twocomp.py --natoms 4170 --skip 1
"""

import argparse
import cPickle
import copy

import numpy as np

class Counts(object):
    """
    Class to take care of the counting for either an acceptor or donor
    Distinguish between trehalose-water (T-W), trehalose-trehalose (T-T)
    and water-water (W-W) hydogren bonds
    """
    def __init__(self, natoms):
        self.ww = 0
        self.tt = 0
        self.tw = 0
        self.total = 0
        self.nsnap = 0
        self.natoms = natoms
        self.curr = {}

    def __add__(self, other) :
        both = Counts(self.natoms)
        both.ww = self.ww + other.ww
        both.tt = self.tt + other.tt
        both.tw = self.tw + other.tw
        both.total = self.total + other.total
        both.nsnap = self.nsnap
        return both

    def __radd__(self, other) :
        both = Counts(self.natoms)
        both.ww = self.ww + other.ww
        both.tt = self.tt + other.tt
        both.tw = self.tw + other.tw
        both.total = self.total + other.total
        both.nsnap = self.nsnap
        return both

    def add_row(self, atomidx, row):
        """
        This will be called when a new bond is found in the
        database. It will check that only chemical hydrogen bonds
        are acounted for, i.e. it will remove double counting
        """
        if atomidx in self.curr:
            #self.curr[-atomidx] = self.curr[atomidx]
            #self.curr[atomidx] = row
            if self.curr[atomidx].distance > row.distance :
                self.curr[atomidx] = row
        else:
            self.curr[atomidx] = row

    def update(self):
        """
        This is to be called every new timestep and classifies the bonds
        """
        self.nsnap += 1
        for id, row in self.curr.iteritems():
            self._classify(row)
        self.curr = {}

    def _classify(self, row):
        """
        This classifies a hydrogen bond as either T-W, T-T or W-W
        """
        if row is None : return
        self.total += 1
        if (row.acceptor_resnm in ["0GA","1GA"] and row.donor_resnm == "SOL") or \
            (row.donor_resnm in ["0GA","1GA"] and row.acceptor_resnm == "SOL") :
            self.tw += 1
        elif row.acceptor_resnm in ["0GA","1GA"] and row.donor_resnm in ["0GA","1GA"] :
            self.tt += 1
        elif row.acceptor_resnm == "SOL"  and row.donor_resnm == "SOL" :
            self.ww += 1
        else:
            print row.acceptor_resnm , row.donor_resnm

    def _reduce(self, lst):
        if len(lst) == 0 :
            return None
        elif len(lst) == 1 :
            return lst[0]
        else :
            idx = np.argmin([row.distance for row in lst])
            return lst[idx]

    def __str__(self):
        fac = 1.0  / float(self.nsnap)
        #return "\n".join([
        #        "\tT-W:\t%.4f"%(self.tw*fac),
        #        "\tT-T:\t%.4f"%(self.tt*fac),
        #        "\tW-W:\t%.4f"%(self.ww*fac),
        #        "\tTotal:\t%.4f"%(self.total*fac),
        #    ])
        return "\t%9.2f\t%9.2f\t%9.2f\t%9.2f"%(self.tw*fac,self.tt*fac,self.ww*fac,self.total*fac)

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description="Script analyse hbond db")
    argparser.add_argument('-f', '--file', help="the pickle file", default="hbonds.pickle")
    argparser.add_argument('-n', '--natoms', type=int, help="the total number of atoms")
    argparser.add_argument('-s', '--skip', type=int, help="the number of snapshots to skip")
    args = argparser.parse_args()

    natoms = args.natoms
    p = cPickle.load(open(args.file))

    # Setup count objects for each donor and acceptor atom
    donors = "H2O H3O H4O H6O HW1 HW2".split()
    donors_counts = {d : Counts(natoms) for d in donors}
    acceptors = "O2 O3 O4 O6 OW".split()
    acceptors_counts = {a : Counts(natoms) for a in acceptors}
    lookup = np.zeros((natoms,natoms),dtype=int)
    id = 0
    for i in range(natoms):
        for j in range(natoms):
            lookup[i,j] = id
            id += 1
    pairs = []
    for d in donors :
        for a in acceptors:
            pairs.append(d+"-"+a)
    pairs_counts = {p : Counts(natoms*natoms)  for p in pairs}

    prevtime = None
    nproc = 0
    print "Starting to analyse each row in the database..."
    for row in p:
        if prevtime is not None and prevtime != row.time :
            nproc += 1
            if nproc > args.skip :
                for d in donors:
                    donors_counts[d].update()
                for a in acceptors:
                    acceptors_counts[a].update()
                for p in pairs:
                    pairs_counts[p].update()
        #if row.donor_resid == row.acceptor_resid :
        #    print row
        donors_counts[row.donor_atom].add_row(row.donor_idx-1, row)
        acceptors_counts[row.acceptor_atom].add_row(row.acceptor_idx-1, row)
        id = lookup[row.donor_idx-1,row.acceptor_idx-1]
        pairs_counts[row.donor_atom+"-"+row.acceptor_atom].curr[id] = row
        prevtime = row.time

    for d in donors:
        donors_counts[d].update()
    for a in acceptors:
        acceptors_counts[a].update()
    for p in pairs:
        pairs_counts[p].update()

    print "Analysed %d timestep \n"%(donors_counts[donors[0]].nsnap)
    print "\n%10s%9s%9s%9s%9s"%("","T-W","T-T","W-W","Tot.")
    total_counts = donors_counts[donors[0]]
    print "-- Donors: --"
    print "%10s"%donors[0],
    print donors_counts[donors[0]]
    for d in donors[1:]:
        total_counts += donors_counts[d]
        print "%10s"%d,
        print donors_counts[d]
    print "%10s"%"Total",
    print total_counts

    print "-- Acceptors: --"
    total_counts = acceptors_counts[acceptors[0]]
    print "%10s"%acceptors[0],
    print acceptors_counts[acceptors[0]]
    for a in acceptors[1:]:
        total_counts += acceptors_counts[a]
        print "%10s"%a,
        print acceptors_counts[a]
    print "%10s"%"Total",
    print total_counts

    print "-- Pairs: --"
    total_counts = pairs_counts[pairs[0]]
    print "%10s"%pairs[0],
    print pairs_counts[pairs[0]]
    for p in pairs[1:]:
        total_counts += pairs_counts[p]
        print "%10s"%p,
        print pairs_counts[p]
    print "%10s"%"Total",
    print total_counts
