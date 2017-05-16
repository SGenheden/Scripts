# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to analyse the RDF of a solute with the centre of liposome

Examples:
  md_liposome_rdf.py -s ref.gro -f sim.xtc --lipids "resname DPPG or DPPC" --solute "resname 5al"
"""

import numpy as np
import MDAnalysis as md
import MDAnalysis.lib.distances as mddist

from sgenlib import moldyn
from sgenlib import mdactions
from sgenlib import pbc

class LiposomeRdfAnalysis(mdactions.TrajectoryAction) :

    def add_arguments(self, parser):
        parser.add_argument('--solute', nargs="+", help="the selectiom mask for solute")
        parser.add_argument('--lipids', help="the selection mask for all lipids")
        parser.add_argument('-o','--out',help="the output",default="lipordf.txt")

    def setup(self, args):
        self.solute = [self.processor.universe.select_atoms(s) for s in args.solute]
        print "Solute selection contains %d atoms in %d residues"%(len(self.solute[0]), len(self.solute[0].residues))
        self.lipids = self.processor.universe.select_atoms(args.lipids)
        print "Lipids selection contains %d residues"%len(self.lipids.residues)
        self.out = args.out
        self.resolution = 1.0
        self.edges = np.arange(0.0, self.processor.universe.dimensions[:3].max()+ \
            self.resolution,self.resolution)
        self.density = np.zeros(self.edges.shape[0]+1)

    def process(self):

        def _accumulate_density(solutecoms, lipidcom) :
            r = np.sqrt(((solutecoms-lipidcom)**2).sum(axis=1))
            for rdig in np.digitize(r, self.edges) :
                self.density[rdig] += 1

        lipidcom = np.median(self.lipids.get_positions() ,axis=0)
        if len(self.solute) == 1 :
            solutecoms = np.asarray([r.center_of_geometry() for r in self.solute.residues])
            _accumulate_density(solutecoms, lipidcom)
        else:
            for sel in self.solute:
                _accumulate_density(sel.get_positions(), lipidcom)


    def finalize(self):
        with open(self.out, "w") as f:
            for d,e in zip(self.density[:-1], self.edges):
                if len(self.solute) == 1 :
                    n = len(self.solute.residues)
                else :
                    n = len(self.solute)*len(self.solute[0])
                gr = float(d) / float(self.processor.nprocessed*n)
                f.write("%.2f %d %.6f\n"%(e,d, gr))

if __name__ == '__main__' :

    processor = moldyn.TrajectoryProcessor("Calculate liposome RDf of a solute")
    analysis = LiposomeRdfAnalysis(processor)
    processor.setup(printargs=True)
    processor.process()
