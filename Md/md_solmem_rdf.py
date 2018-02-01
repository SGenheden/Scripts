# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to calculate the RDF of a solute in a membrane

Experimental code!
"""

import numpy as np
import MDAnalysis.lib.distances as mddist
import scipy.spatial.distance as scidist

from sgenlib import moldyn
from sgenlib import mdactions
from sgenlib import pbc

class SolMemRdf(mdactions.TrajectoryAction) :

    def add_arguments(self, parser):
        parser.add_argument('--ref',help="the selectiom mask for the reference")
        parser.add_argument('--sol', nargs="+", help="the selectiom mask for the solute")
        parser.add_argument('-o','--out',help="the output",default="rdf.txt")
        parser.add_argument('--zero',choices=["None","xy"],help="if to zero some components",default="xy")
        parser.add_argument('--resolution', type=float, help="the resolution of the RDF", default=0.5)
        parser.add_argument('--expand', action="store_true", help="expand the RDF around zero", default=False)
        parser.add_argument('--peratom', action="store_true", help="calculate per atom densities", default=False)

    def setup(self, args):
        self.reference = self.processor.universe.select_atoms(args.ref)
        print "Reference selection contains %d atoms in %d residues"%(len(self.reference), len(self.reference.residues))
        self.solute = [self.processor.universe.select_atoms(s) for s in args.sol]
        print "Solute selection contains %d atoms in %d residues"%(len(self.solute[0]), len(self.solute[0].residues))
        self.out = args.out
        self.peratom = args.peratom
        self.expand = args.expand
        self.zero = args.zero
        # Setup density
        self.resolution = args.resolution
        self.edges = np.arange(0.0, self.processor.universe.dimensions[:3].max()*0.5+ \
            self.resolution,self.resolution)
        self.density = np.zeros(self.edges.shape[0]+1)

    def process(self):

        def _accumulate_density(solutecoms, refcom) :
            r = np.sqrt(((solutecoms - refcom)**2).sum(axis=1))
            for rdig in np.digitize(r, self.edges) :
                self.density[rdig] += 1

        refcom = np.median(self.reference.get_positions() , axis=0)
        if self.zero == "xy":
            refcom[0:2] = 0.0
        if len(self.solute) == 1 :
            if not self.peratom :
                solutecoms = np.asarray([r.center_of_geometry() for r in self.solute[0].residues])
                if self.zero == "xy":
                    solutecoms[:, 0:2] = 0.0
                _accumulate_density(solutecoms, refcom)
            else :
                p = self.solute[0].get_positions()
                if self.zero == "xy":
                    p[:, 0:2] = 0.0
                _accumulate_density(p, refcom)
        else:
            for sel in self.solute:
                p = sel.get_positions()
                if self.zero == "xy":
                    p[:, 0:2] = 0.0
                _accumulate_density(p, refcom)

    def finalize(self):
        if len(self.solute) == 1 :
            if not self.peratom:
                n = len(self.solute[0].residues)
            else:
                n = len(self.solute[0])
        else :
            n = len(self.solute)*len(self.solute[0])
        norm = 1 / float(self.processor.nprocessed*n*self.resolution)
        with open(self.out, "w") as f:
            if self.expand :
                for d, e, in zip(self.density[::-1][1:], -self.edges[::-1]):
                    f.write("%.2f %d %.6f\n"%(e, d, float(d) * norm))
            for d, e in zip(self.density[:-1], self.edges):
                f.write("%.2f %d %.6f\n"%(e, d, float(d) * norm))

if __name__ == '__main__' :

    processor = moldyn.TrajectoryProcessor("Calculate the RDF of solute in a membrane")
    analysis = SolMemRdf(processor)
    processor.setup(printargs=True)
    processor.process()
