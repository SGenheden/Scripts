# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to analyze the extent of alcohols
"""

import numpy as np
import MDAnalysis.core.distances as mddist

from sgenlib import moldyn
from sgenlib import mdactions
from sgenlib import pbc

class AlcoholLenAnalysis(mdactions.TrajectoryAction) :

    def add_arguments(self, parser):
        parser.add_argument('--sel',nargs=2,help="the selectiom mask for the atoms")
        parser.add_argument('-o','--out',help="the output",default="alcohol_len.txt")

    def setup(self, args):
        self.sel1 = self.processor.universe.select_atoms(args.sel[0])
        self.sel2 = self.processor.universe.select_atoms(args.sel[1])
        print "Group selection contains %d and %d atoms"%(len(self.sel1), len(self.sel2))
        self.out = args.out
        self.records = []

    def process(self):

        zlen = np.mean(np.abs(self.sel1.positions[:, 2] -
                                self.sel2.positions[:, 2]))
        diff2 = np.power(self.sel1.positions - self.sel2.positions, 2)
        totlen = np.mean(np.sqrt(np.sum(diff2, axis=1)))
        self.records.append(mdactions.MDRecord(self.processor.currtime,[totlen, zlen]))

    def finalize(self):
        self._write_records(headers="Tot_len Z_len".split())

if __name__ == '__main__' :

    processor = moldyn.TrajectoryProcessor("Measure length of alcohol")
    analysis = AlcoholLenAnalysis(processor)
    processor.setup(printargs=True)
    processor.process()
