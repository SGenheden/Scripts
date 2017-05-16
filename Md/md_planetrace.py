# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to analyse the trace of atoms in a plane

Examples:
  md_planetrace.py -f md.xtc -s md.tpr --sel "resname AAC and name C1"
"""

import numpy as np
import matplotlib.pylab as plt

from sgenlib import moldyn
from sgenlib import mdactions
from sgenlib import colors

class PlaneTraceAnalysis(mdactions.TrajectoryAction) :

    def add_arguments(self, parser):
        parser.add_argument('--sel',help="the selectiom mask")
        parser.add_argument('-o','--out',help="the output",default="planetrace")

    def setup(self, args):
        self.atoms = self.processor.universe.select_atoms(args.sel)
        self.box = self.processor.universe.dimensions[:2]*0.1
        self.halfbox = self.box / 2.0
        print "%d atoms selected"%len(self.atoms)
        self.traces = [[] for atom in self.atoms]
        self.outpre = args.out
        self.records = []

    def process(self):
        for i, pos in enumerate(self.atoms.positions):
            self.traces[i].append(mdactions.MDRecord(
                    self.processor.currtime, pos[:2]*0.1))

    def finalize(self):
        ncols = 2
        nrows = int(np.ceil(len(self.traces)/float(ncols)))
        f = plt.figure()
        for i, trace in enumerate(self.traces, 1):
            a = f.add_subplot(nrows, ncols, i)
            self._plot_trace(trace, a, i-1)

            self.records = trace
            self.out = self.outpre + "-%d.txt"%i
            self._write_records()

        f.savefig(self.outpre + ".png", format="png")

    def _plot_trace(self, trace, axis, idx):
        """
        Plot the trace, but make breaks if the
        atoms pass over the periodic box so that very long lines aren't drawn
        """
        xy = [trace[0].value]
        for prevrecord, record in zip(trace[:-1], trace[1:]):
            if np.abs(prevrecord.value[0]-record.value[0]) > self.halfbox[0] or \
                np.abs(prevrecord.value[1]-record.value[1]) > self.halfbox[1] :
                xy = np.asarray(xy)
                axis.plot(xy[:,0], xy[:,1], "-", color=colors.color(idx))
                xy = [record.value]
            else:
                xy.append(record.value)
        xy = np.asarray(xy)
        axis.plot(xy[:,0], xy[:,1], "-", color=colors.color(idx))
        axis.set_xlim(0, self.box[0])
        axis.set_ylim(0, self.box[1])
        axis.set_aspect('equal')

if __name__ == '__main__' :

    processor = moldyn.TrajectoryProcessor("Calculate the plane trace of atoms")
    analysis = PlaneTraceAnalysis(processor)
    processor.setup(printargs=True)
    processor.process()
