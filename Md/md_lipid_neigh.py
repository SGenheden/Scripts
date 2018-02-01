# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to analyse the neighbors of lipids

Examples:

"""

import numpy as np
import MDAnalysis as md
import MDAnalysis.lib.distances as mddist
import scipy.spatial.distance as scidist

from sgenlib import moldyn
from sgenlib import mdactions
from sgenlib import pbc

class LipidNeighborAnalysis(mdactions.TrajectoryAction) :

    def add_arguments(self, parser):
        parser.add_argument('--head', help="the selection mask for the heads")
        parser.add_argument('--lipid', help="the selection mask for all lipids")
        parser.add_argument('--nneigh', type=int, help="the number of neightbors to consider", default=10)
        parser.add_argument('-o','--out',help="the output",default="lipidneigh.txt")

    def setup(self, args):
        self.lipid = self.processor.universe.select_atoms(args.lipid)
        print "Lipids selection contains %d residues"%len(self.lipid.residues)
        self.head = self.processor.universe.select_atoms(args.head)
        print "Head selection contains %d particles"%len(self.head)
        self.out = args.out
        self.nneigh = args.nneigh
        self.resolution = 0.5
        self.samples = []

    def process(self):

        def _add_distance(positions, box):
            dist = scidist.squareform(mddist.self_distance_array(positions, box))
            for i, d in enumerate(dist):
                self.samples.extend(np.sort(d)[1:self.nneigh+1])

        com = np.asarray([r.center_of_geometry() for r in self.lipid.residues])
        midz = self.head.positions[:,2].mean()
        lowsel = com[:,2] < midz
        uppsel = np.logical_not(lowsel)
        com[:,2] = 0.0
        _add_distance(com[lowsel,:], self.processor.currsnap.dimensions)
        _add_distance(com[uppsel,:], self.processor.currsnap.dimensions)


    def finalize(self):
        self.samples = np.asarray(self.samples)
        edges = np.arange(self.samples.min(), self.samples.max()+ \
            self.resolution, self.resolution)
        his, b = np.histogram(self.samples, edges)
        dens, b = np.histogram(self.samples, edges, density=True)
        mid = 0.5 * (edges[1:]+edges[:-1])
        with open(self.out, "w") as f:
            f.write("# Mean = %.3f \n# Std = %.3f \n# Median = %.3f\n"%(self.samples.mean(),self.samples.std(),np.median(self.samples)))
            for e, h, d  in zip(mid, his, dens):
                f.write("%.2f %d %.6f\n"%(e, h, d))

if __name__ == '__main__' :

    processor = moldyn.TrajectoryProcessor("Calculate neighbor distance distribution")
    analysis = LipidNeighborAnalysis(processor)
    processor.setup(printargs=True)
    processor.process()
